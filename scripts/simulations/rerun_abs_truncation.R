#!/usr/bin/env Rscript
# ============================================================================
# Re-run of the simulations whose joint-4PL ABSOLUTE CTmax_1hr was mis-read
# (2026-11-26).
#
# What went wrong: the retired bayesTLS::extract_tdt() read the absolute
# (literal 50% survival) CTmax_1hr by inverting each draw's LT curve on a fixed
# grid, min(T) - 2 to max(T) + 2 C. A draw whose crossing lay outside the grid
# was dropped from the summary without a warning, and a curve that is not
# monotone over the grid could get a temperature that is not a crossing at all.
# extract_joint() (sim_functions.R) fixes both; nothing else changes. Draws
# that still have no single falling crossing within 30 C (for example an LT
# that lengthens with temperature) stay out of the summary; they are counted
# per simulation (meta$abs_uncrossed) and per cell in the before/after table.
#
# How: the June software (a frozen bayesTLS that has the fix, CmdStan 2.38.0)
# reproduces the joint fits exactly, so each simulation is re-run from its own
# seed. The relative arm and the absolute z must come out identical. The
# two-stage rows are kept from the June files (merge_rerun()): that arm has no
# bug, and the glmmTMB installed now (1.1.15.2) differs from June's on a few
# datasets. Steps: back up the raw files -> re-run into a staging folder ->
# merge and check -> write back -> re-aggregate -> before/after summary.
# Resume-safe: a simulation already re-run (staging + merged file) is skipped.
# ============================================================================
PKG_LIB <- file.path(Sys.getenv("HOME"), "Library", "Caches", "tls-sim-lib")  # PROVENANCE.txt inside
CMDSTAN <- file.path(Sys.getenv("HOME"), ".cmdstan", "cmdstan-2.38.0")       # the June runs' CmdStan
Sys.setenv(CMDSTAN = CMDSTAN)          # read by cmdstanr when it loads (workers inherit it)
.libPaths(c(PKG_LIB, .libPaths()))
suppressPackageStartupMessages({
  library(bayesTLS); library(dplyr); library(tibble); library(tidyr); library(here)
})
stopifnot(cmdstanr::cmdstan_version() == "2.38.0",
          normalizePath(find.package("bayesTLS")) == normalizePath(file.path(PKG_LIB, "bayesTLS")))
# Master and workers source one frozen copy of sim_functions.R, made on the
# first run, so editing the repo copy mid-run cannot change the cells to come.
FROZEN_SRC <- here::here("output", "sim_twostage", "rerun_stage_2026-11-26",
                         "sim_functions_frozen.R")
if (!file.exists(FROZEN_SRC)) {
  dir.create(dirname(FROZEN_SRC), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(here::here("scripts", "simulations", "sim_functions.R"), FROZEN_SRC))
}
source(FROZEN_SRC)
source(here::here("scripts", "simulations", "sim_scenarios.R"))

# ---- 1. Settings ------------------------------------------------------------
OUT_DIR     <- here::here("output", "sim_twostage")
RAW_ROOT    <- file.path(OUT_DIR, "raw")
BACKUP_ROOT <- file.path(OUT_DIR, "raw_pre_rerun_2026-11-26")   # the June files, untouched
STAGE_ROOT  <- file.path(OUT_DIR, "rerun_stage_2026-11-26")     # re-run results + checks
MASTER_SEED <- 20260513L                                        # as run_simulations.R
SAMPLER     <- list(chains = 3, iter = 3000, warmup = 1500, max_treedepth = 16, adapt_delta = 0.95)
WORKERS     <- 4L
# Which simulations: in check_rerun_scope.R (first 6 sims of every cell) the
# absolute CTmax changed only in simulations that had also lost draws, the
# non-crossing values included. So every simulation that lost a draw is
# re-run; a cell where at least HEAVY of its simulations did is re-run in full;
# in the other cells N_VERIFY randomly chosen simulations that lost no draw are
# re-run as a check, and if any of them changes the whole cell is re-run.
HEAVY    <- 0.10
N_VERIFY <- 30L
N_SIMS   <- 1000L                      # simulations per cell, as run_simulations.R
dir.create(STAGE_ROOT, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Which simulations ------------------------------------------------------
scope <- readRDS(here::here("output", "sim_twostage_check", "rerun_scope.rds"))
lost <- lapply(setNames(scenarios$label, scenarios$label), function(lab) {
  f <- list.files(file.path(RAW_ROOT, lab), "^sim_\\d+\\.rds$", full.names = TRUE)
  ids <- as.integer(sub("^sim_(\\d+)\\.rds$", "\\1", basename(f)))
  # Read the June copy where a backup exists: a restart then sees the files as
  # they were before the re-run overwrote them, so the plan is rebuilt exactly.
  bak <- file.path(BACKUP_ROOT, lab, basename(f))
  f <- ifelse(file.exists(bak), bak, f)
  ids[vapply(f, function(x) { r <- readRDS(x); isTRUE(r$meta$joint_ok) && NROW(r$draws_abs) < 1000 },
             logical(1))]
})
scope_changed_clean <- scope |> dplyr::filter(abs_ct_d > 1e-8, kept_old == 1000)
CELLS <- union(scenarios$label[lengths(lost) > 0], unique(scope_changed_clean$scenario))
set.seed(20261126L)                    # the verification sample
plan <- lapply(setNames(CELLS, CELLS), function(lab) {
  all_ids <- seq_len(N_SIMS)
  heavy <- length(lost[[lab]]) >= HEAVY * N_SIMS || lab %in% scope_changed_clean$scenario
  clean <- setdiff(all_ids, lost[[lab]])
  list(heavy = heavy, lost = lost[[lab]],
       verify = if (heavy) integer(0) else sort(sample(clean, min(N_VERIFY, length(clean)))),
       ids = if (heavy) all_ids else sort(c(lost[[lab]], NULL)))
})
for (lab in CELLS) plan[[lab]]$ids <- sort(union(plan[[lab]]$ids, plan[[lab]]$verify))
print(data.frame(cell = CELLS, sims_lost_draws = vapply(plan, function(x) length(x$lost), 1L),
                 rerun_in_full = vapply(plan, `[[`, NA, "heavy"),
                 sims_rerun = vapply(plan, function(x) length(x$ids), 1L)), row.names = FALSE)
saveRDS(plan, file.path(STAGE_ROOT, "rerun_plan.rds"))

# ---- 3. Back up, re-run, merge, write back -----------------------------------
for (lab in CELLS) {
  sc <- scenarios[scenarios$label == lab, ]
  truth <- sim_truth(dgp = sc$dgp, family = sc$family, design = sc$design, u_0 = na_to_null(sc$u_0),
                     ell_0 = na_to_null(sc$ell_0), u_beta1 = na_to_null(sc$u_beta1))
  raw_dir <- file.path(RAW_ROOT, lab); bak_dir <- file.path(BACKUP_ROOT, lab)
  stage_dir <- file.path(STAGE_ROOT, lab)
  dir.create(bak_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(raw_dir, "^sim_\\d+\\.rds$")
  todo_bak <- setdiff(files, list.files(bak_dir, "^sim_\\d+\\.rds$"))
  if (length(todo_bak))
    stopifnot(all(file.copy(file.path(raw_dir, todo_bak), file.path(bak_dir, todo_bak))))
  stopifnot(length(list.files(bak_dir, "^sim_\\d+\\.rds$")) == length(files))

  run_one <- function(sim_id) {   # re-run one sim, keep the June two-stage rows, check
    f <- sprintf("sim_%05d.rds", sim_id)
    if (file.exists(file.path(stage_dir, sub("sim_", "check_", f)))) return(invisible(NULL))
    run_one_sim(sim_id, sc, truth, stage_dir, SAMPLER, MASTER_SEED, force = TRUE)
    new <- readRDS(file.path(stage_dir, f))
    if (!isTRUE(new$meta$joint_ok)) return(invisible(NULL))  # failed fit: June file kept, no check
    mg <- merge_rerun(new, readRDS(file.path(bak_dir, f)))
    save_retry(mg$result, file.path(raw_dir, f))
    save_retry(mg$check, file.path(stage_dir, sub("sim_", "check_", f)))
    invisible(NULL)
  }
  run_ids <- function(ids) {
    todo <- ids[!file.exists(file.path(stage_dir, sprintf("check_%05d.rds", ids)))]
    message(sprintf("%s: %d of %d sims to re-run", lab, length(todo), length(ids)))
    if (length(todo))
      run_parallel(todo, run_one, workers = WORKERS, envir = environment(), lib = PKG_LIB,
                   src = FROZEN_SRC,
                   exports = c("sc", "truth", "raw_dir", "bak_dir", "stage_dir", "SAMPLER",
                               "MASTER_SEED", "run_one"))
    left <- ids[!file.exists(file.path(stage_dir, sprintf("check_%05d.rds", ids)))]
    if (length(left))
      stop(lab, ": no successful re-run for sims ", toString(left),
           "; run the script again to retry them", call. = FALSE)
  }
  run_ids(plan[[lab]]$ids)
  # Verification: if any re-run simulation that lost no draw changed, re-run the cell.
  if (!plan[[lab]]$heavy) {
    chk <- dplyr::bind_rows(lapply(file.path(stage_dir, sprintf("check_%05d.rds", plan[[lab]]$verify)), readRDS))
    n_chg <- sum(is.na(chk$abs_ct_d) | chk$abs_ct_d > 1e-8)
    message(sprintf("%s: verification sample, %d of %d simulations without lost draws changed",
                    lab, n_chg, nrow(chk)))
    if (n_chg > 0L) { plan[[lab]]$heavy <- TRUE; run_ids(seq_len(N_SIMS)) }
  }
}
saveRDS(plan, file.path(STAGE_ROOT, "rerun_plan.rds"))

# ---- 3b. Merge every re-run simulation again ----------------------------------
# raw = merge_rerun(re-run, June backup) for every simulation with a check, so a
# change to merge_rerun() after a run reaches the files that run wrote. A
# re-run whose fit failed (an earlier version merged those) gets its June file
# back and loses its check, so the next run retries it.
for (lab in CELLS) {
  stage_dir <- file.path(STAGE_ROOT, lab); raw_dir <- file.path(RAW_ROOT, lab)
  bak_dir   <- file.path(BACKUP_ROOT, lab)
  ids <- as.integer(sub("^check_(\\d+)\\.rds$", "\\1",
                        list.files(stage_dir, "^check_\\d+\\.rds$")))
  failed <- integer(0)
  for (id in ids) {
    f   <- sprintf("sim_%05d.rds", id)
    new <- readRDS(file.path(stage_dir, f))
    if (!isTRUE(new$meta$joint_ok)) {
      stopifnot(file.copy(file.path(bak_dir, f), file.path(raw_dir, f), overwrite = TRUE))
      unlink(file.path(stage_dir, sub("sim_", "check_", f)))
      failed <- c(failed, id); next
    }
    mg <- merge_rerun(new, readRDS(file.path(bak_dir, f)))
    save_retry(mg$result, file.path(raw_dir, f))
    save_retry(mg$check, file.path(stage_dir, sub("sim_", "check_", f)))
  }
  if (length(failed))
    stop(lab, ": re-run fits failed for sims ", toString(failed),
         "; their June files are restored; run the script again to retry them", call. = FALSE)
}

# ---- 4. Checks: everything but the absolute CTmax must reproduce --------------
checks <- dplyr::bind_rows(lapply(CELLS, function(lab) {
  f <- list.files(file.path(STAGE_ROOT, lab), "^check_\\d+\\.rds$", full.names = TRUE)
  dplyr::bind_rows(lapply(f, readRDS))
}))
saveRDS(checks, file.path(STAGE_ROOT, "rerun_checks.rds"))
bad <- dplyr::filter(checks, is.na(rel_maxd) | is.na(abs_z_d) | rel_maxd > 1e-6 | abs_z_d > 1e-6)
message(sprintf("re-run sims: %d; relative arm or absolute z NOT reproduced in %d", nrow(checks), nrow(bad)))
if (nrow(bad)) {
  print(as.data.frame(bad))
  stop("the re-run did not reproduce; cells not re-aggregated", call. = FALSE)
}

# ---- 5. Re-aggregate the re-run cells ----------------------------------------
for (lab in CELLS)
  print(as.data.frame(aggregate_scenario(lab, file.path(RAW_ROOT, lab), OUT_DIR)))

# ---- 6. Before / after, absolute CTmax_1hr -----------------------------------
before_after <- dplyr::bind_rows(lapply(CELLS, function(lab) {
  old <- collect_raw(file.path(BACKUP_ROOT, lab))$per_sim
  new <- readRDS(file.path(OUT_DIR, sprintf("per_sim_%s.rds", lab)))
  sm <- function(d) d |> dplyr::filter(method == "joint_4pl_abs", quantity == "CTmax_1hr") |>
    # rmse before bias: summarise() evaluates in order, so a later `bias^2`
    # would read the summarised mean, not the per-simulation errors
    dplyr::summarise(coverage = mean(covered, na.rm = TRUE), rmse = sqrt(mean(bias^2, na.rm = TRUE)),
                     bias = mean(bias, na.rm = TRUE), width = mean(width, na.rm = TRUE))
  dplyr::bind_cols(tibble::tibble(scenario = lab), dplyr::rename_with(sm(old), ~ paste0(.x, "_before")),
                   dplyr::rename_with(sm(new), ~ paste0(.x, "_after")))
}))
# Per cell, among the re-run simulations: June fits that dropped draws, and
# fits that still leave draws without a single falling crossing (left out of
# their summary; see the header).
unc <- checks |>
  dplyr::group_by(scenario) |>
  dplyr::summarise(sims_dropped_june   = sum(kept_old < 1000),
                   sims_uncrossed_now  = sum(uncrossed > 0, na.rm = TRUE),
                   draws_uncrossed_now = sum(uncrossed, na.rm = TRUE),
                   max_uncrossed_now   = max(uncrossed, na.rm = TRUE), .groups = "drop")
before_after <- dplyr::left_join(before_after, unc, by = "scenario")
saveRDS(before_after, file.path(OUT_DIR, "rerun_2026-11-26_before_after.rds"))
# A copy in output/data (the OSF results tier) for the supplement's table of
# draws left out at the absolute threshold.
saveRDS(before_after, here::here("output", "data", "sim_rerun_before_after.rds"))
print(as.data.frame(before_after), digits = 4)
