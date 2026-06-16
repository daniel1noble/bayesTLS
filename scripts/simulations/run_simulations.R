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

suppressPackageStartupMessages({
  library(bayesTLS); library(dplyr); library(tibble)
  library(tidyr);    library(parallel); library(here)
})
source(here::here("scripts", "simulations", "sim_functions.R"))

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

# ---- 2. Scenarios — the single source of truth (replaces the bash queue) ---
# Each row is plain data: a label plus the sim_truth() arguments for that cell.
# NA in u_0 / ell_0 / u_beta1 means "keep the DGP preset". Edit a number to
# change a sweep value; add a row to add a cell.
#
#   scen1  strict-equivalence baseline (binomial DGM, asymptotes at the bound)
#   scen2  likelihood misspecification only (beta-binomial, same shape)
#   scen3  heat lowers max survival      (asymmetric u drift)
#   scen4  asymptotes compress           (symmetric u + ell drift)
#   scen5  heat sharpens the curve       (k varies with T)
#   scen6  sweep: strength of u drift    (u_beta1)
#   scen7  sweep: level of upper asymptote (u_0)
#   scen8  sweep: design intensity       (full/sparse x replication)
#   scen9  sweep: exposure-time window   (max duration 60/120/240/405 min) -- the
#          design question: how well does a shortened assay recover z & CTmax_1hr?
scenarios <- tibble::tribble(
  ~label,                   ~dgp,        ~family,         ~n_reps, ~design,   ~u_0,  ~ell_0, ~u_beta1,
  "scen1_strict_eq_n3",     "baseline",  "binomial",      3,  "full",    0.999, 0.001,  NA,
  "scen1_strict_eq_n5",     "baseline",  "binomial",      5,  "full",    0.999, 0.001,  NA,
  "scen2_lik_misspec_n3",   "baseline",  "beta_binomial", 3,  "full",    0.999, 0.001,  NA,
  "scen2_lik_misspec_n5",   "baseline",  "beta_binomial", 5,  "full",    0.999, 0.001,  NA,
  "scen3_heat_lowers_u_n3", "asym_u",    "beta_binomial", 3,  "full",    NA,    NA,     NA,
  "scen3_heat_lowers_u_n5", "asym_u",    "beta_binomial", 5,  "full",    NA,    NA,     NA,
  "scen4_compress_n3",      "sym_ul",    "beta_binomial", 3,  "full",    NA,    NA,     NA,
  "scen4_compress_n5",      "sym_ul",    "beta_binomial", 5,  "full",    NA,    NA,     NA,
  "scen5_sharpen_n3",       "varying_k", "beta_binomial", 3,  "full",    NA,    NA,     NA,
  "scen5_sharpen_n5",       "varying_k", "beta_binomial", 5,  "full",    NA,    NA,     NA,
  "scen6_ub_m005",          "baseline",  "beta_binomial", 5,  "full",    NA,    NA,     -0.005,
  "scen6_ub_m010",          "baseline",  "beta_binomial", 5,  "full",    NA,    NA,     -0.010,
  "scen6_ub_m015",          "baseline",  "beta_binomial", 5,  "full",    NA,    NA,     -0.015,
  "scen6_ub_m019",          "baseline",  "beta_binomial", 5,  "full",    NA,    NA,     -0.019,
  "scen7_u0_099",           "baseline",  "beta_binomial", 5,  "full",    0.99,  NA,     NA,
  "scen7_u0_095",           "baseline",  "beta_binomial", 5,  "full",    0.95,  NA,     NA,
  "scen7_u0_085",           "baseline",  "beta_binomial", 5,  "full",    0.85,  NA,     NA,
  "scen7_u0_075",           "baseline",  "beta_binomial", 5,  "full",    0.75,  NA,     NA,
  "scen7_u0_065",           "baseline",  "beta_binomial", 5,  "full",    0.65,  NA,     NA,
  "scen8_full_n1",          "baseline",  "beta_binomial", 1,  "full",    NA,    NA,     NA,
  "scen8_full_n3",          "baseline",  "beta_binomial", 3,  "full",    NA,    NA,     NA,
  "scen8_full_n5",          "baseline",  "beta_binomial", 5,  "full",    NA,    NA,     NA,
  "scen8_sparse_n1",        "baseline",  "beta_binomial", 1,  "sparse",  NA,    NA,     NA,
  "scen8_sparse_n3",        "baseline",  "beta_binomial", 3,  "sparse",  NA,    NA,     NA,
  "scen8_sparse_n5",        "baseline",  "beta_binomial", 5,  "sparse",  NA,    NA,     NA,
  "scen9_tmax_060",         "baseline",  "beta_binomial", 5,  "tmax060", NA,    NA,     NA,
  "scen9_tmax_120",         "baseline",  "beta_binomial", 5,  "tmax120", NA,    NA,     NA,
  "scen9_tmax_240",         "baseline",  "beta_binomial", 5,  "tmax240", NA,    NA,     NA,
  "scen9_tmax_405",         "baseline",  "beta_binomial", 5,  "tmax405", NA,    NA,     NA,
)
scenarios$index <- seq_len(nrow(scenarios))   # disjoint seed offset per scenario

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

  # The whole per-simulation story, in plain sight:
  #   simulate one dataset -> fit two-stage (both Stage-1 likelihoods) ->
  #   fit the joint 4PL -> score all three against the truth -> save.
  # Resume-safe: a sim whose raw file exists is skipped unless FORCE.
  run_one <- function(sim_id) {
    f <- file.path(raw_dir, sprintf("sim_%05d.rds", sim_id))
    if (!FORCE && file.exists(f)) return(invisible(NULL))
    seed   <- MASTER_SEED + 1000L * sc$index + sim_id
    data   <- sim_dataset(n_reps = sc$n_reps, seed = seed, truth = truth)
    ts_bin <- fit_two_stage(data, stage1 = "binomial")
    ts_bb  <- fit_two_stage(data, stage1 = "betabinomial")
    joint  <- do.call(fit_joint_4pl, c(list(data, seed = seed), SAMPLER))
    save_retry(score_run(joint, ts_bin, ts_bb, truth, sim_id, sc$label, seed), f)
    invisible(NULL)
  }

  todo <- setdiff(seq_len(N_SIMS), done_ids(raw_dir, FORCE))
  message(sprintf("  %d sims to run (%d cached)", length(todo), N_SIMS - length(todo)))
  if (length(todo))
    run_parallel(todo, run_one, workers = WORKERS,
                 exports = c("sc", "truth", "raw_dir", "SAMPLER",
                             "MASTER_SEED", "FORCE", "run_one"))

  # Aggregate the raw per-sim files into the 8 scenario objects. Each save is
  # one named helper call, so the output set is trivial to track.
  agg <- collect_raw(raw_dir)
  out <- function(prefix, obj)
    save_retry(obj, file.path(OUT_DIR, sprintf("%s_%s.rds", prefix, sc$label)))
  out("per_sim",      agg$per_sim)
  out("meta",         agg$meta)
  out("draws",        agg$draws)
  out("draws_abs",    agg$draws_abs)
  out("summary",      summarise_mcse(agg$per_sim))
  out("summary_conv", summarise_mcse(agg$per_sim, conv = TRUE, meta = agg$meta))
  out("diag",         diag_summary(agg$meta, sc$label))
  out("diffs",        pairwise_diffs(agg$per_sim))

  print(as.data.frame(summarise_mcse(agg$per_sim)))
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
