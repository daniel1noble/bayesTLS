#!/usr/bin/env Rscript
# Two-stage bias simulation driver.
#
# Generates N_SIMS beta-binomial 4PL datasets, fits both the joint Bayesian
# 4PL and the classical two-stage TDT pipeline to each, and saves each sim's
# result to its own RDS file. Resume-safe: re-running skips any sim whose
# result file already exists (unless --force is set).
#
# Parallelism uses a PSOCK cluster (`parallel::makeCluster` +
# `parallel::parLapply`). Fresh R workers; we explicitly export the few
# script-level objects each worker needs and load `bayesTLS` once per worker.
# A PSOCK cluster is needed (instead of `mclapply`'s fork) because cmdstanr
# launches Stan model executables via processx, and processx is not
# fork-safe on macOS. Set --workers = 1 for serial debugging.
#
# Truth: beta-binomial 4PL with temperature only on `mid`. Truth values and
# the (T × duration) design come from bayesTLS::sim_twostage_truth() and
# bayesTLS::sim_twostage_grid().
#
# CLI:
#   Rscript scripts/sim_twostage_bias.R --scenario n3 --n_sims 1000 \
#       --workers 5 --out_dir output/sim_twostage --seed 20260513
#
# Run both scenarios with two launches.

suppressPackageStartupMessages({
  library(bayesTLS)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(parallel)
  library(optparse)
  library(here)
})

# Simulation-specific helpers live alongside this driver, not inside the
# bayesTLS package. They depend on bayesTLS's public API (standardize_data,
# fit_4pl, extract_tdt, diagnose_tdt_fit) but are private to this manuscript's
# bias-simulation study.
source(here::here("scripts", "sim_twostage_helpers.R"))

# ---- CLI ---------------------------------------------------------------------

option_list <- list(
  optparse::make_option("--scenario", type = "character", default = "n3",
                        help = "Replication-budget shorthand: n3 or n5 (sets n_reps to 3 or 5). Ignored when --n_reps is given."),
  optparse::make_option("--n_reps",   type = "integer",  default = NA_integer_,
                        help = "Replicate cups per (T, t) cell. Overrides --scenario."),
  optparse::make_option("--dgp",      type = "character", default = "baseline",
                        help = "DGP preset: baseline, sym_ul, asym_u, or varying_k."),
  optparse::make_option("--design",   type = "character", default = "full",
                        help = "Design grid: 'full' (5T x 6t) or 'sparse' (3T x 4t)."),
  optparse::make_option("--u_0",       type = "double", default = NA_real_,
                        help = "Override upper asymptote at T_bar. Default = DGP preset."),
  optparse::make_option("--u_beta1",   type = "double", default = NA_real_,
                        help = "Override temperature slope of u. Default = DGP preset."),
  optparse::make_option("--ell_beta1", type = "double", default = NA_real_,
                        help = "Override temperature slope of ell. Default = DGP preset."),
  optparse::make_option("--k_beta1",   type = "double", default = NA_real_,
                        help = "Override temperature slope of k. Default = DGP preset."),
  optparse::make_option("--label",    type = "character", default = NA_character_,
                        help = "Custom output label (file/dir name). Default: derived from --dgp/--scenario."),
  optparse::make_option("--n_sims",   type = "integer",  default = 1000L,
                        help = "Number of simulated datasets."),
  optparse::make_option("--workers",  type = "integer",  default = 5L,
                        help = "Parallel workers (PSOCK cluster size)."),
  optparse::make_option("--out_dir",  type = "character",
                        default = "output/sim_twostage",
                        help = "Output directory for per-sim and summary RDS."),
  optparse::make_option("--seed",     type = "integer",  default = 20260513L,
                        help = "Master seed; sim seed = master + 1000*scen_id + sim_id."),
  optparse::make_option("--force",    type = "logical",  default = FALSE,
                        help = "If TRUE, re-run sims even if result file exists.")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

stopifnot(opt$dgp %in% c("baseline", "sym_ul", "asym_u", "varying_k"))
stopifnot(opt$design %in% c("full", "sparse"))

# Resolve n_reps. --n_reps overrides --scenario; otherwise --scenario must be
# n3/n5. Both routes set scen_id (used to keep simulation seeds disjoint).
if (!is.na(opt$n_reps)) {
  n_reps  <- opt$n_reps
  scen_id <- as.integer(100L + n_reps)
} else {
  stopifnot(opt$scenario %in% c("n3", "n5"))
  n_reps_map <- list(n3 = 3L, n5 = 5L)
  n_reps     <- n_reps_map[[opt$scenario]]
  scen_id    <- match(opt$scenario, names(n_reps_map))
}

master <- opt$seed

# Convert NA overrides → NULL for sim_twostage_truth(), which uses NULL to
# mean "keep the DGP preset value".
to_null <- function(x) if (is.na(x)) NULL else x
truth  <- sim_twostage_truth(
  dgp       = opt$dgp,
  u_0       = to_null(opt$u_0),
  u_beta1   = to_null(opt$u_beta1),
  ell_beta1 = to_null(opt$ell_beta1),
  k_beta1   = to_null(opt$k_beta1),
  design    = opt$design
)

# Output label: explicit --label wins; else built from the DGP + scenario.
# (opt$scenario always has a default of "n3", so the baseline branch is
# self-contained.)
path_label <- if (!is.na(opt$label)) {
  opt$label
} else if (opt$dgp == "baseline") {
  opt$scenario
} else {
  paste(opt$dgp, opt$scenario, sep = "_")
}
stopifnot(is.character(path_label), nchar(path_label) > 0)

raw_dir <- file.path(opt$out_dir, "raw", path_label)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("\n=== Two-stage bias simulation ===\n"))
cat(sprintf("  label:       %s\n", path_label))
cat(sprintf("  dgp:         %s\n", opt$dgp))
cat(sprintf("  design:      %s\n", opt$design))
cat(sprintf("  n_reps:      %d (scenario = %s)\n", n_reps, opt$scenario))
cat(sprintf("  n_sims:      %d\n", opt$n_sims))
cat(sprintf("  workers:     %d\n", opt$workers))
cat(sprintf("  out_dir:     %s\n", opt$out_dir))
cat(sprintf("  raw_dir:     %s\n", raw_dir))
cat(sprintf("  master seed: %d\n", master))
cat(sprintf("  force:       %s\n\n", opt$force))
cat(sprintf("  truth: u_0=%g  ell_0=%g  k_0=%g\n",
            truth$u, truth$ell, truth$k))
cat(sprintf("  truth slopes: u_beta1=%g  ell_beta1=%g  k_beta1=%g\n",
            truth$u_beta1, truth$ell_beta1, truth$k_beta1))
cat(sprintf("  OLS-derived: z_true = %.4f °C, CTmax_1hr_true = %.4f °C\n\n",
            truth$z_true, truth$CTmax_1hr_true))

# ---- one simulation ----------------------------------------------------------
# Forked workers inherit `n_reps`, `scen_id`, `master`, `truth`, `raw_dir`,
# and the `opt$force` flag from the parent — nothing to wire up.

run_one <- function(sim_id) {
  result_file <- file.path(raw_dir, sprintf("sim_%05d.rds", sim_id))
  if (!opt$force && file.exists(result_file)) return(invisible(NULL))

  seed_sim <- master + 1000L * scen_id + sim_id
  data <- sim_twostage_dataset(n_reps = n_reps, seed = seed_sim)

  t0 <- Sys.time()
  joint <- tryCatch(
    fit_joint_4pl_sim(data, seed = seed_sim),
    error = function(e) list(success = FALSE, error = conditionMessage(e),
                              z = list(point = NA_real_, lower = NA_real_, upper = NA_real_),
                              CTmax_1hr = list(point = NA_real_, lower = NA_real_,
                                               upper = NA_real_))
  )
  joint_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ts <- fit_two_stage_classical(data)

  row  <- sim_twostage_result_row(joint, ts, truth, sim_id, path_label,
                                  runtime_sec = joint_sec)
  meta <- tibble::tibble(
    sim_id      = sim_id,
    scenario    = path_label,
    seed        = seed_sim,
    joint_sec   = joint_sec,
    joint_ok    = joint$success,
    ts_ok       = ts$success,
    rhat_max    = if (joint$success) joint$diagnostics$rhat_max else NA_real_,
    divergences = if (joint$success) joint$diagnostics$divergences else NA_integer_
  )

  # Full per-draw joint-4PL posterior (z, CTmax_1hr, T_crit). Tagged with
  # sim_id so downstream code can pool draws across sims without losing
  # the per-sim grouping. The two-stage point estimates and Wald CIs are
  # already in `row` — no per-draw analogue for the classical pipeline.
  draws_df <- if (joint$success && !is.null(joint$draws)) {
    joint$draws |>
      dplyr::mutate(sim_id   = sim_id,
                    scenario = path_label) |>
      dplyr::select(sim_id, scenario, .draw, z, CTmax_1hr, T_crit)
  } else NULL

  saveRDS(list(row = row, meta = meta, draws = draws_df), result_file)
  invisible(NULL)
}

# ---- run (resume-safe) -------------------------------------------------------

existing <- list.files(raw_dir, pattern = "^sim_\\d+\\.rds$")
done_ids <- if (length(existing) == 0L || opt$force) integer(0) else
  as.integer(sub("^sim_(\\d+)\\.rds$", "\\1", existing))
to_run   <- setdiff(seq_len(opt$n_sims), done_ids)
cat(sprintf("  %d sims to run (%d already complete)\n\n",
            length(to_run), length(done_ids)))

t_start <- Sys.time()
if (length(to_run) > 0L) {
  if (opt$workers > 1L) {
    cl <- parallel::makeCluster(opt$workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    helpers_path <- here::here("scripts", "sim_twostage_helpers.R")
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages({
      library(bayesTLS); library(dplyr); library(tibble)
    }))
    parallel::clusterCall(cl, source, helpers_path)
    parallel::clusterExport(
      cl,
      varlist = c("opt", "n_reps", "scen_id", "master", "truth",
                   "raw_dir", "path_label", "run_one"),
      envir = environment()
    )
    invisible(parallel::parLapply(cl, to_run, run_one))
  } else {
    invisible(lapply(to_run, run_one))
  }
}
elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat(sprintf("\n  Run elapsed: %.1f min (%.2f sec / sim avg with %d workers)\n",
            elapsed / 60,
            if (length(to_run) > 0L) elapsed / length(to_run) else 0,
            opt$workers))

# ---- aggregate per-sim files -------------------------------------------------

cat("\n  Aggregating per-sim files ...\n")
all_files <- list.files(raw_dir, pattern = "^sim_\\d+\\.rds$", full.names = TRUE)
all_rows  <- vector("list", length(all_files))
all_meta  <- vector("list", length(all_files))
all_draws <- vector("list", length(all_files))
for (k in seq_along(all_files)) {
  obj <- readRDS(all_files[k])
  all_rows[[k]]  <- obj$row
  all_meta[[k]]  <- obj$meta
  all_draws[[k]] <- obj$draws
}
per_sim <- dplyr::bind_rows(all_rows)
meta    <- dplyr::bind_rows(all_meta)
draws   <- dplyr::bind_rows(all_draws)

if (nrow(per_sim) == 0L) {
  stop("No per-sim files were produced — every simulation failed. ",
       "Check the worker error messages above (typically Stan/cmdstanr issues).")
}

# ---- summary with Monte Carlo standard errors --------------------------------

mcse_summary <- per_sim |>
  dplyr::filter(success) |>
  dplyr::group_by(scenario, method, quantity) |>
  dplyr::summarise(
    n            = dplyr::n(),
    mean_bias    = mean(bias, na.rm = TRUE),
    mcse_bias    = stats::sd(bias, na.rm = TRUE) / sqrt(dplyr::n()),
    rmse         = sqrt(mean(bias^2, na.rm = TRUE)),
    mcse_rmse    = stats::sd(bias^2, na.rm = TRUE) /
                    (2 * sqrt(mean(bias^2, na.rm = TRUE)) * sqrt(dplyr::n())),
    coverage     = mean(covered, na.rm = TRUE),
    mcse_cov     = sqrt(mean(covered, na.rm = TRUE) *
                         (1 - mean(covered, na.rm = TRUE)) / dplyr::n()),
    med_width    = stats::median(width, na.rm = TRUE),
    .groups      = "drop"
  ) |>
  dplyr::mutate(dplyr::across(c(mean_bias, mcse_bias, rmse, mcse_rmse,
                                 coverage, mcse_cov, med_width),
                              ~ round(.x, 4)))

# ---- paired method-difference summary (joint_4pl - two_stage per sim) -------
# These come straight from the per-sim table by pivoting on `method` and
# subtracting. Captures whether the two methods agree on each individual
# dataset (different from comparing aggregate biases, which marginalise out
# the per-sim correlation).

# Paired diff: keep only sims where BOTH methods succeeded (otherwise the
# pivot is missing a column and a per-sim contrast is undefined). When the
# data are very sparse (e.g. sparse design x n_reps=1) all two-stage fits
# may fail; in that case `diffs` stays empty and we record a placeholder.
both_ok <- per_sim |>
  dplyr::select(sim_id, scenario, method, quantity, estimate, success) |>
  tidyr::pivot_wider(names_from  = method,
                     values_from = c(estimate, success),
                     names_sep   = "_")

needed_cols <- c("estimate_joint_4pl", "estimate_two_stage",
                 "success_joint_4pl",  "success_two_stage")
diffs <- if (all(needed_cols %in% names(both_ok))) {
  both_ok |>
    dplyr::filter(success_joint_4pl, success_two_stage) |>
    dplyr::transmute(sim_id, scenario, quantity,
                     joint_4pl = estimate_joint_4pl,
                     two_stage = estimate_two_stage,
                     diff      = joint_4pl - two_stage)
} else {
  tibble::tibble(sim_id = integer(0), scenario = character(0),
                 quantity = character(0),
                 joint_4pl = numeric(0), two_stage = numeric(0),
                 diff = numeric(0))
}

diff_summary <- if (nrow(diffs) == 0L) {
  tibble::tibble(scenario = character(0), quantity = character(0),
                 n = integer(0), mean_diff = numeric(0),
                 mcse_diff = numeric(0), median_diff = numeric(0),
                 diff_q025 = numeric(0), diff_q975 = numeric(0))
} else {
  diffs |>
    dplyr::group_by(scenario, quantity) |>
    dplyr::summarise(
      n            = dplyr::n(),
      mean_diff    = mean(diff, na.rm = TRUE),
      mcse_diff    = stats::sd(diff, na.rm = TRUE) / sqrt(dplyr::n()),
      median_diff  = stats::median(diff, na.rm = TRUE),
      diff_q025    = stats::quantile(diff, 0.025, na.rm = TRUE, names = FALSE),
      diff_q975    = stats::quantile(diff, 0.975, na.rm = TRUE, names = FALSE),
      .groups      = "drop"
    ) |>
    dplyr::mutate(dplyr::across(c(mean_diff, mcse_diff, median_diff,
                                   diff_q025, diff_q975),
                                ~ round(.x, 4)))
}

out_per_sim <- file.path(opt$out_dir, sprintf("per_sim_%s.rds", path_label))
out_meta    <- file.path(opt$out_dir, sprintf("meta_%s.rds",    path_label))
out_draws   <- file.path(opt$out_dir, sprintf("draws_%s.rds",   path_label))
out_summary <- file.path(opt$out_dir, sprintf("summary_%s.rds", path_label))
out_diffs   <- file.path(opt$out_dir, sprintf("diffs_%s.rds",   path_label))
saveRDS(per_sim,      out_per_sim)
saveRDS(meta,         out_meta)
saveRDS(draws,        out_draws)
saveRDS(mcse_summary, out_summary)
saveRDS(diffs,        out_diffs)
cat(sprintf("\n  Saved %s (%d rows)\n", out_per_sim, nrow(per_sim)))
cat(sprintf("  Saved %s (%d rows)\n", out_meta, nrow(meta)))
cat(sprintf("  Saved %s (%d rows of per-draw joint-4PL posterior)\n",
            out_draws, nrow(draws)))
cat(sprintf("  Saved %s\n", out_summary))
cat(sprintf("  Saved %s (per-sim joint-vs-two-stage differences)\n", out_diffs))

cat("\n=== Joint-4PL minus Two-stage paired difference ===\n")
print(as.data.frame(diff_summary))

cat("\n=== Summary with Monte Carlo SEs ===\n")
print(as.data.frame(mcse_summary))

# ---- MCSE decision rule ------------------------------------------------------

cat("\n=== MCSE decision check ===\n")
cat("  Thresholds: coverage MCSE < 0.015, bias MCSE < 0.05 °C\n\n")
check <- mcse_summary |>
  dplyr::mutate(
    cov_ok  = mcse_cov  < 0.015,
    bias_ok = mcse_bias < 0.05
  )
print(as.data.frame(check[, c("scenario", "method", "quantity",
                               "mcse_bias", "bias_ok",
                               "mcse_cov", "cov_ok")]))

all_pass <- all(check$cov_ok) && all(check$bias_ok)
if (all_pass) {
  cat(sprintf("\n  All MCSEs below thresholds at n_sims = %d. Stopping here.\n",
              opt$n_sims))
} else {
  cat(sprintf("\n  Some MCSEs above thresholds at n_sims = %d.\n", opt$n_sims))
  cat("  Bump --n_sims by 500 and re-run; existing per-sim files will be reused.\n")
}

cat("\nDONE.\n")
