#!/usr/bin/env Rscript
# ============================================================================
# Scope check for the 2026-11-26 re-run (see rerun_abs_truncation.R).
#
# Re-runs the first N_CHECK simulations of EVERY scenario with the fixed joint
# extraction (extract_joint()) on the June software (frozen bayesTLS + CmdStan
# 2.38.0) into a scratch folder, and compares each with its stored raw file.
# The relative arm and the absolute z must reproduce exactly (they are
# untouched); the absolute CTmax_1hr changes where the retired extract_tdt()
# dropped draws off its grid, or returned a non-crossing for a curve that is
# not monotone over the grid. Cells where nothing changes can be left alone.
# Writes output/sim_twostage_check/rerun_scope.rds (gitignored).
# ============================================================================
PKG_LIB <- file.path(Sys.getenv("HOME"), "Library", "Caches", "tls-sim-lib")
CMDSTAN <- file.path(Sys.getenv("HOME"), ".cmdstan", "cmdstan-2.38.0")   # the June runs' CmdStan
Sys.setenv(CMDSTAN = CMDSTAN)          # read by cmdstanr when it loads (workers inherit it)
.libPaths(c(PKG_LIB, .libPaths()))
suppressPackageStartupMessages({
  library(bayesTLS); library(dplyr); library(tibble); library(tidyr); library(here)
})
stopifnot(cmdstanr::cmdstan_version() == "2.38.0")
source(here::here("scripts", "simulations", "sim_functions.R"))
source(here::here("scripts", "simulations", "sim_scenarios.R"))

N_CHECK     <- 6L
WORKERS     <- 3L
MASTER_SEED <- 20260513L
SAMPLER     <- list(chains = 3, iter = 3000, warmup = 1500, max_treedepth = 16, adapt_delta = 0.95)
RAW_ROOT    <- here::here("output", "sim_twostage", "raw")
SCRATCH     <- here::here("output", "sim_twostage_check", "rerun_scope")

compare_sim <- function(new, old) {
  key <- c("method", "quantity")
  j <- merge(old$row[, c(key, "estimate", "lower", "upper")], new$row[, c(key, "estimate", "lower", "upper")],
             by = key, suffixes = c(".old", ".new"))
  d <- pmax(abs(j$estimate.new - j$estimate.old), abs(j$lower.new - j$lower.old), abs(j$upper.new - j$upper.old))
  pick <- function(m, q) d[j$method == m & j$quantity == q]
  tibble::tibble(rel_maxd   = max(d[j$method == "joint_4pl"]),
                 abs_z_d    = pick("joint_4pl_abs", "z"),
                 abs_ct_d   = pick("joint_4pl_abs", "CTmax_1hr"),
                 two_stage_maxd = max(d[grepl("^two_stage", j$method)]),
                 kept_old = nrow(old$draws_abs), kept_new = nrow(new$draws_abs),
                 uncrossed = new$meta$abs_uncrossed)
}

out <- list()
for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]
  truth <- sim_truth(dgp = sc$dgp, family = sc$family, design = sc$design, u_0 = na_to_null(sc$u_0),
                     ell_0 = na_to_null(sc$ell_0), u_beta1 = na_to_null(sc$u_beta1))
  raw_dir <- file.path(SCRATCH, sc$label); dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  run_one <- function(sim_id) run_one_sim(sim_id, sc, truth, raw_dir, SAMPLER, MASTER_SEED, force = TRUE)
  run_parallel(seq_len(N_CHECK), run_one, workers = WORKERS, envir = environment(), lib = PKG_LIB,
               exports = c("sc", "truth", "raw_dir", "SAMPLER", "MASTER_SEED", "run_one"))
  out[[sc$label]] <- dplyr::bind_rows(lapply(seq_len(N_CHECK), function(id) {
    f <- sprintf("sim_%05d.rds", id)
    dplyr::bind_cols(tibble::tibble(scenario = sc$label, sim_id = id),
                     compare_sim(readRDS(file.path(raw_dir, f)), readRDS(file.path(RAW_ROOT, sc$label, f))))
  }))
  message(sprintf("%-24s rel max|d| %.1e  abs z %.1e  abs CTmax changed in %d of %d sims",
                  sc$label, max(out[[sc$label]]$rel_maxd), max(out[[sc$label]]$abs_z_d),
                  sum(out[[sc$label]]$abs_ct_d > 1e-8), N_CHECK))
}
scope <- dplyr::bind_rows(out)
saveRDS(scope, file.path(dirname(SCRATCH), "rerun_scope.rds"))
print(as.data.frame(scope), digits = 3)
