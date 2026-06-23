#!/usr/bin/env Rscript
# ============================================================================
# Honest re-validation of the simulation against the reworked bayesTLS package.
#
# Re-runs the first N simulations of selected scenarios with the CURRENT package
# and compares them to the pre-rework snapshot (output/sim_twostage_snapshot_*).
# The new functions run exactly as written -- NOTHING is tuned to reproduce the
# old numbers. Two layers of comparison:
#
#   * two_stage_*  -- deterministic glm/glmmTMB, independent of bayesTLS. Because
#     sim_id fixes the seed (and sim_dataset/sim_truth are unchanged), the same
#     data goes in, so these MUST match the snapshot exactly. They are a built-in
#     CONTROL: any drift here means the reconstruction is off, not the package.
#   * joint_4pl*   -- MCMC fits through fit_4pl()+extract_tdt() (the entire change
#     surface). Compared on per-sim point-estimate drift AND operating
#     characteristics (coverage), since sampling jitter means individual sims are
#     never bit-identical even with no code change. A near-zero MEAN signed drift
#     = sampling jitter; a systematic non-zero mean = a real change to chase down.
#
# Writes ONLY to output/sim_twostage_check/ (scratch). Never touches the snapshot
# or output/sim_twostage/.
#
#   Rscript scripts/simulations/check_against_snapshot.R                       # defaults
#   Rscript scripts/simulations/check_against_snapshot.R 25                    # N = 25
#   Rscript scripts/simulations/check_against_snapshot.R 25 scen1_strict_eq_n3 # one scenario
# ============================================================================
suppressPackageStartupMessages({
  library(bayesTLS); library(dplyr); library(tibble); library(here); library(parallel)
})
source(here::here("scripts", "simulations", "sim_functions.R"))

SNAP_DIR    <- here::here("output", "sim_twostage_snapshot_2026-06-23")
SCRATCH     <- here::here("output", "sim_twostage_check")
MASTER_SEED <- 20260513L                       # must match run_simulations.R
WORKERS     <- 5L
SAMPLER     <- list(chains = 3, iter = 3000, warmup = 1500,
                    max_treedepth = 16, adapt_delta = 0.95)

# Scenario rows mirror run_simulations.R EXACTLY, including `index` (the row
# number there) which sets the per-scenario seed offset. Do not change these.
SCEN <- tibble::tribble(
  ~label,                   ~index, ~dgp,       ~family,         ~n_reps, ~design,  ~u_0,  ~ell_0, ~u_beta1,
  "scen1_strict_eq_n3",     1L,     "baseline",  "binomial",      3, "full",   0.999, 0.001,  NA,
  "scen3_heat_lowers_u_n3", 5L,     "asym_u",    "beta_binomial", 3, "full",   NA,    NA,     NA,
  "scen5_sharpen_n3",       9L,     "varying_k", "beta_binomial", 3, "full",   NA,    NA,     NA,
  "scen8_sparse_n3",        24L,    "baseline",  "beta_binomial", 3, "sparse", NA,    NA,     NA
)

# Robustness: the `index` (seed offset) MUST equal each scenario's row number in
# run_simulations.R, or the reconstructed seeds desync and we compare different
# datasets. The two_stage control catches a desync after the fact, but verify it
# up front by deriving the canonical order straight from the runner and asserting
# the hand-typed indices agree.
.runner_labels <- local({
  ln <- readLines(here::here("scripts", "simulations", "run_simulations.R"))
  m  <- regmatches(ln, regexpr('"scen[0-9][a-z0-9_]*"', ln))   # quoted tribble labels, in order
  gsub('"', '', m[nzchar(m)])
})
.derived <- match(SCEN$label, .runner_labels)
if (anyNA(.derived) || !identical(as.integer(.derived), as.integer(SCEN$index)))
  stop("SCEN$index disagrees with run_simulations.R order:\n",
       paste(sprintf("  %-24s manual=%s derived=%s", SCEN$label, SCEN$index, .derived),
             collapse = "\n"), call. = FALSE)

# ---- args ------------------------------------------------------------------
args   <- commandArgs(trailingOnly = TRUE)
N      <- if (length(args) && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 16L
labels <- args[grepl("^scen", args)]
if (!length(labels)) labels <- c("scen1_strict_eq_n3", "scen3_heat_lowers_u_n3")
SCEN   <- dplyr::filter(SCEN, label %in% labels)
stopifnot(nrow(SCEN) > 0)

# ---- run one scenario's first N sims into the scratch dir ------------------
run_scenario_check <- function(sc) {
  truth <- sim_truth(dgp = sc$dgp, family = sc$family, design = sc$design,
                     u_0 = na_to_null(sc$u_0), ell_0 = na_to_null(sc$ell_0),
                     u_beta1 = na_to_null(sc$u_beta1))
  raw_dir <- file.path(SCRATCH, "raw", sc$label)
  unlink(raw_dir, recursive = TRUE)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  run_one <- function(sim_id) {
    f      <- file.path(raw_dir, sprintf("sim_%05d.rds", sim_id))
    seed   <- MASTER_SEED + 1000L * sc$index + sim_id      # identical to the runner
    data   <- sim_dataset(n_reps = sc$n_reps, seed = seed, truth = truth)
    ts_bin <- fit_two_stage(data, stage1 = "binomial")
    ts_bb  <- fit_two_stage(data, stage1 = "betabinomial")
    joint  <- do.call(fit_joint_4pl, c(list(data, seed = seed), SAMPLER))
    save_retry(score_run(joint, ts_bin, ts_bb, truth, sim_id, sc$label, seed), f)
    invisible(NULL)
  }
  message(sprintf("[%s] running %d sims (index=%d, z_true=%.3f, CTmax_true=%.3f)...",
                  sc$label, N, sc$index, truth$z_true, truth$CTmax_1hr_true))
  run_parallel(seq_len(N), run_one, workers = WORKERS,
               exports = c("sc", "truth", "raw_dir", "SAMPLER", "MASTER_SEED", "run_one"),
               envir = environment())
  collect_raw(raw_dir)$per_sim
}

# ---- run + compare ---------------------------------------------------------
dir.create(SCRATCH, recursive = TRUE, showWarnings = FALSE)
cmp_rows <- list()
for (i in seq_len(nrow(SCEN))) {
  sc  <- SCEN[i, ]
  new <- run_scenario_check(sc)
  old <- readRDS(file.path(SNAP_DIR, sprintf("per_sim_%s.rds", sc$label))) |>
    dplyr::filter(sim_id <= N)
  j <- dplyr::inner_join(
    dplyr::transmute(new, sim_id, method, quantity,
                     est_new = estimate, lo_new = lower, hi_new = upper, cov_new = covered),
    dplyr::transmute(old, sim_id, method, quantity,
                     est_old = estimate, lo_old = lower, hi_old = upper, cov_old = covered),
    by = c("sim_id", "method", "quantity"))
  cmp_rows[[sc$label]] <- j |>
    dplyr::group_by(method, quantity) |>
    dplyr::summarise(
      n              = dplyr::n(),
      mean_signed_d  = mean(est_new - est_old, na.rm = TRUE),   # systematic shift?
      max_abs_d_est  = max(abs(est_new - est_old), na.rm = TRUE),
      max_abs_d_ci   = max(abs(lo_new - lo_old), abs(hi_new - hi_old), na.rm = TRUE),
      cov_old        = mean(cov_old, na.rm = TRUE),
      cov_new        = mean(cov_new, na.rm = TRUE),
      .groups = "drop") |>
    dplyr::mutate(scenario = sc$label, .before = 1)
}
cmp <- dplyr::bind_rows(cmp_rows)

cat("\n================ OLD (snapshot) vs NEW (reworked pkg) ================\n")
cat("CONTROL: two_stage_* must show max_abs_d_est ~0 (deterministic).\n")
cat("SIGNAL : joint_4pl* -- mean_signed_d ~0 => MCMC jitter; non-zero => real shift.\n\n")
print(as.data.frame(cmp), digits = 4)
saveRDS(list(comparison = cmp, N = N, scenarios = SCEN$label),
        file.path(SCRATCH, "comparison_summary.rds"))
cat("\nsaved:", file.path(SCRATCH, "comparison_summary.rds"), "\n")
