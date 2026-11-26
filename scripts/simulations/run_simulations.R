#!/usr/bin/env Rscript
# ============================================================================
# Two-stage-bias simulation — single runner.
#
# Everything the simulation does is visible here or in the building blocks it
# sources from sim_functions.R. No shell dispatcher, no CLI parser. To follow
# one simulation, read run_one() below: simulate -> two-stage x2 -> joint 4PL
# -> score. To follow the whole study, read the `scenarios` table.
#
# Run ALL scenarios:
#   Rscript scripts/simulations/run_simulations.R
# Run ONE (or a few) by label:
#   Rscript scripts/simulations/run_simulations.R scen9_tmax_060 scen9_tmax_405
# Or set SCENARIOS_TO_RUN below and source() the file interactively.
# ============================================================================

# Optional private library holding the bayesTLS version to run (NULL = the
# default library). The 2026-11-26 re-run used a frozen copy (see
# rerun_abs_truncation.R) so that package edits during the run could not reach it.
PKG_LIB <- NULL
if (!is.null(PKG_LIB)) .libPaths(c(PKG_LIB, .libPaths()))

# CmdStan version. The June 2026 runs used 2.38.0, and reproducing their joint
# fits needs it (2.40.0 gives different draws). Set before cmdstanr loads; NULL
# = cmdstanr's default installation (then the MCMC draws differ from June's).
# Two-stage rows also depend on the glmmTMB version: see PROVENANCE.md.
CMDSTAN_VERSION <- "2.38.0"
if (!is.null(CMDSTAN_VERSION)) {
  cmdstan_dir <- file.path(Sys.getenv("HOME"), ".cmdstan", paste0("cmdstan-", CMDSTAN_VERSION))
  if (!dir.exists(cmdstan_dir))
    stop("CmdStan ", CMDSTAN_VERSION, " not found at ", cmdstan_dir, ". Install it with ",
         "cmdstanr::install_cmdstan(version = \"", CMDSTAN_VERSION, "\"), or set ",
         "CMDSTAN_VERSION <- NULL to use the default installation.", call. = FALSE)
  Sys.setenv(CMDSTAN = cmdstan_dir)
}

suppressPackageStartupMessages({
  library(bayesTLS); library(dplyr); library(tibble)
  library(tidyr);    library(parallel); library(here)
})
source(here::here("scripts", "simulations", "sim_functions.R"))
if (!is.null(CMDSTAN_VERSION)) stopifnot(cmdstanr::cmdstan_version() == CMDSTAN_VERSION)

# ---- 1. Config (edit these) ------------------------------------------------
OUT_DIR     <- here::here("output", "sim_twostage")
# Raw per-sim .rds files live in-repo under OUT_DIR/raw (gitignored) so the whole
# workflow runs from the project root: clone the repo, drop the raw files from
# OSF into output/sim_twostage/raw/, then re-aggregate or re-run. The durable
# copy of the ~29k raw files lives on OSF, not Git. NOTE: on a Dropbox-synced
# checkout the file provider may dehydrate these many small files to online-only;
# rehydrate them (Finder > Make Available Offline) or re-pull from OSF before
# aggregation. Point RAW_ROOT at a non-synced path to opt out of Dropbox.
RAW_ROOT    <- file.path(OUT_DIR, "raw")
N_SIMS      <- 1000          # simulated datasets per scenario
WORKERS     <- 5             # PSOCK cluster size (cmdstanr is not fork-safe)
MASTER_SEED <- 20260513L     # seed_sim = MASTER_SEED + 1000*scenario_index + sim_id
FORCE       <- FALSE         # TRUE = ignore the raw/ cache and refit every sim

# Sampler settings handed straight to fit_joint_4pl().
SAMPLER <- list(chains = 3, iter = 3000, warmup = 1500,
                max_treedepth = 16, adapt_delta = 0.95)

# NULL = run every scenario; trailing CLI args (or an inline vector) subset it.
cli <- commandArgs(trailingOnly = TRUE)
SCENARIOS_TO_RUN <- if (length(cli) > 0) cli else NULL

# ---- 2. Scenarios -----------------------------------------------------------
# The scenario table lives in sim_scenarios.R (shared with the re-run script);
# its row order sets each scenario's seed offset.
source(here::here("scripts", "simulations", "sim_scenarios.R"))

# ---- 3. One scenario, start to finish --------------------------------------
# Builds the truth, runs N_SIMS datasets through all three estimators, writes
# the 8 output objects (same filenames the manuscript figures read).
run_scenario <- function(sc) {
  truth <- sim_truth(dgp = sc$dgp, family = sc$family, design = sc$design,
                     u_0 = na_to_null(sc$u_0), ell_0 = na_to_null(sc$ell_0),
                     u_beta1 = na_to_null(sc$u_beta1))
  message(sprintf("\n=== %s  (z_true=%.3f  CTmax_1hr_true=%.3f) ===",
                  sc$label, truth$z_true, truth$CTmax_1hr_true))

  raw_dir <- file.path(RAW_ROOT, sc$label)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  # The whole per-simulation story is run_one_sim() in sim_functions.R:
  #   simulate one dataset -> fit two-stage (both Stage-1 likelihoods) ->
  #   fit the joint 4PL -> score all three against the truth -> save.
  # Resume-safe: a sim whose raw file exists is skipped unless FORCE.
  run_one <- function(sim_id)
    run_one_sim(sim_id, sc, truth, raw_dir, SAMPLER, MASTER_SEED, FORCE)

  todo <- setdiff(seq_len(N_SIMS), done_ids(raw_dir, FORCE))
  message(sprintf("  %d sims to run (%d cached)", length(todo), N_SIMS - length(todo)))
  if (length(todo))
    run_parallel(todo, run_one, workers = WORKERS,
                 exports = c("sc", "truth", "raw_dir", "SAMPLER",
                             "MASTER_SEED", "FORCE", "run_one"), lib = PKG_LIB)

  # Aggregate the raw per-sim files into the 8 scenario objects the documents
  # read (aggregate_scenario() in sim_functions.R names each one).
  print(as.data.frame(aggregate_scenario(sc$label, raw_dir, OUT_DIR)))
  invisible(NULL)
}

# ---- 4. Pre-flight, then run ------------------------------------------------
to_run <- if (is.null(SCENARIOS_TO_RUN)) {
  scenarios
} else {
  dplyr::filter(scenarios, label %in% SCENARIOS_TO_RUN)
}
if (nrow(to_run) == 0L)
  stop("No scenarios matched: ", paste(SCENARIOS_TO_RUN, collapse = ", "))

message(sprintf("Pre-flighting %d scenario(s) ...", nrow(to_run)))
preflight(to_run)

t0 <- Sys.time()
for (i in seq_len(nrow(to_run))) run_scenario(to_run[i, ])
message(sprintf("\nDONE — %d scenario(s) in %.1f min.",
                nrow(to_run), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
