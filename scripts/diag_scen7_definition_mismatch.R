# Definition-mismatch diagnostic for Scen 7 u0=0.65.
#
# Refits a small batch of joint-4PL fits using the saved seeds from
# output/sim_twostage/raw/scen7_u0_065/, then extracts CTmax_1hr both
# ways (target_surv = "relative" vs "absolute"). Prints per-sim bias
# under each definition against the analytical absolute truth.

suppressPackageStartupMessages({
  library(dplyr)
  library(bayesTLS)
  source(here::here("scripts", "sim_twostage_helpers.R"))
})

n_sims_test <- 20  # Refit count (each ~20-25 s at u0=0.65)
ndraws      <- 1000

# Reconstruct the Scen 7 u0=0.65 truth exactly as the driver did.
truth <- sim_twostage_truth(
  dgp    = "baseline",
  u_0    = 0.65,
  family = "beta_binomial"
)
cat(sprintf("Truth (analytical, absolute): z = %.4f, CTmax_1hr = %.4f\n\n",
            truth$z_true, truth$CTmax_1hr_true))

# Pull the seeds from the first n_sims_test sims (deterministic from saved meta).
raw_dir <- here::here("output", "sim_twostage", "raw", "scen7_u0_065")
sim_files <- sort(list.files(raw_dir, full.names = TRUE))[seq_len(n_sims_test)]
seeds <- vapply(sim_files,
                function(f) readRDS(f)$meta$seed,
                FUN.VALUE = integer(1))

results <- vector("list", length(seeds))
for (i in seq_along(seeds)) {
  s <- seeds[i]
  cat(sprintf("[%2d/%d] seed = %d ", i, length(seeds), s))
  t0 <- Sys.time()

  data <- sim_twostage_dataset(n_reps = 5L, seed = s, truth = truth)
  wf   <- standardize_data(data,
                           temp = "T", duration = "t",
                           n_total = "n", n_surv = "y",
                           duration_unit = "minutes")
  wf   <- fit_4pl(wf,
                  chains = 2, iter = 2000, cores = 2,
                  seed = s, refresh = 0, silent = 2,
                  backend = "cmdstanr",
                  control = list(adapt_delta = 0.95, max_treedepth = 14))

  et_rel <- suppressMessages(extract_tdt(
    wf, target_surv = "relative", t_ref = 60,
    time_multiplier = 1, ndraws = ndraws, lethal = TRUE
  ))
  et_abs <- suppressMessages(extract_tdt(
    wf, target_surv = "absolute", t_ref = 60,
    time_multiplier = 1, ndraws = ndraws, lethal = TRUE
  ))

  ctmax_rel <- et_rel$CTmax$summary$temp_median
  ctmax_abs <- et_abs$CTmax$summary$temp_median
  z_rel     <- et_rel$z$summary$z_median
  z_abs     <- et_abs$z$summary$z_median

  results[[i]] <- tibble::tibble(
    sim       = i,
    seed      = s,
    z_rel     = z_rel,
    z_abs     = z_abs,
    ctmax_rel = ctmax_rel,
    ctmax_abs = ctmax_abs,
    bias_rel  = ctmax_rel - truth$CTmax_1hr_true,
    bias_abs  = ctmax_abs - truth$CTmax_1hr_true,
    runtime_s = as.numeric(Sys.time() - t0, units = "secs")
  )
  cat(sprintf("CTmax rel=%.3f abs=%.3f  bias_rel=%+.3f bias_abs=%+.3f  (%.1fs)\n",
              ctmax_rel, ctmax_abs, results[[i]]$bias_rel,
              results[[i]]$bias_abs, results[[i]]$runtime_s))
}

R <- dplyr::bind_rows(results)
out_path <- here::here("output", "sim_twostage",
                       "diag_scen7_u0_065_definition_check.rds")
saveRDS(R, out_path)
cat("\nSaved per-sim results to ", out_path, "\n", sep = "")

cat("\n=== Summary across ", nrow(R), " refit sims ===\n", sep = "")
cat(sprintf("Truth CTmax (absolute): %.4f °C\n", truth$CTmax_1hr_true))
cat(sprintf("Mean CTmax (relative extract): %.4f  →  bias = %+.4f °C\n",
            mean(R$ctmax_rel), mean(R$bias_rel)))
cat(sprintf("Mean CTmax (absolute extract): %.4f  →  bias = %+.4f °C\n",
            mean(R$ctmax_abs), mean(R$bias_abs)))
cat(sprintf("\nMean per-sim offset (absolute - relative): %+.4f °C\n",
            mean(R$ctmax_abs - R$ctmax_rel)))
cat(sprintf("\nFor comparison, the cell's headline bias was +0.6709 °C (n=500, relative extract).\n"))
