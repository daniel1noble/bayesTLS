# One-sim walk-through diagnostic.
#
# Loads ONE per-sim record from output/sim_twostage/raw/<scenario>/sim_<id>.rds,
# regenerates the dataset deterministically from the saved seed, refits all
# three estimators (joint Bayesian 4PL via brms, per-T binomial GLM, per-T
# beta-binomial GLMM), prints each fitted model, computes z and CTmax_1hr
# under both `target_surv = "relative"` and `"absolute"`, then prints a
# side-by-side comparison against the values stored in the original sim row
# to confirm exact reproducibility.
#
# Configure the cell to inspect by editing the constants below. Each refit
# takes ~5–25 s depending on scenario (Scen 1 and Scen 7 cells are slowest).
#
# Per-sim raw .rds layout (from sim_twostage_bias.R:224):
#   $row    one row per (method, quantity): truth, estimate, lower, upper,
#           bias, covered, width, success, runtime_sec
#   $meta   sim_id, scenario, seed, joint_sec, all convergence diagnostics
#   $draws  one row per posterior draw: .draw, z, CTmax_1hr, T_crit
#           (extracted with target_surv = "relative" — the simulation default)

suppressPackageStartupMessages({
  library(dplyr)
  library(bayesTLS)
  library(glmmTMB)
  source(here::here("scripts", "sim_twostage_helpers.R"))
})

# ============================================================================
# Configuration — edit these to walk through a different sim.
# ============================================================================
SCENARIO   <- "scen7_u0_065"
SIM_ID     <- 1L
N_REPS     <- 5L
TRUTH_ARGS <- list(dgp = "baseline", u_0 = 0.65, family = "beta_binomial")

# ============================================================================
# Step 1 — Load saved per-sim record from raw/.
# ============================================================================
raw_path <- here::here("output", "sim_twostage", "raw", SCENARIO,
                       sprintf("sim_%05d.rds", SIM_ID))
saved <- readRDS(raw_path)

cat("\n============================================================\n")
cat(sprintf("Step 1 — Loaded %s\n", raw_path))
cat("============================================================\n")
cat(sprintf("Scenario: %s   Sim id: %d   Saved seed: %d\n",
            saved$meta$scenario, saved$meta$sim_id, saved$meta$seed))
cat("\nSaved per-method point estimates (the numbers the simulation recorded):\n")
print(saved$row[, c("method", "quantity", "truth", "estimate",
                    "lower", "upper", "covered", "success")])

# ============================================================================
# Step 2 — Reconstruct the analytical truth and the dataset (deterministic
#          from the seed). The truth here is "absolute" (50% survival LT50,
#          OLS-projected onto the design temperatures).
# ============================================================================
cat("\n============================================================\n")
cat("Step 2 — Reconstruct truth and dataset from the saved seed\n")
cat("============================================================\n")
truth <- do.call(sim_twostage_truth, TRUTH_ARGS)
cat(sprintf("Truth (analytical, absolute): z = %.4f, CTmax_1hr = %.4f °C\n",
            truth$z_true, truth$CTmax_1hr_true))

data <- sim_twostage_dataset(n_reps = N_REPS, seed = saved$meta$seed,
                             truth = truth)
cat(sprintf("Dataset: %d rows | temperatures = {%s} | durations = {%s}\n",
            nrow(data),
            paste(sort(unique(data$T)), collapse = ", "),
            paste(sort(unique(data$t)), collapse = ", ")))
cat("\nFirst 6 rows of the simulated dataset:\n")
print(head(data, 6))

# ============================================================================
# Step 3 — Refit the joint Bayesian 4PL via brms (cmdstanr backend).
# ============================================================================
cat("\n============================================================\n")
cat("Step 3 — Joint Bayesian 4PL fit (brms, cmdstanr backend)\n")
cat("============================================================\n")
wf <- standardize_data(data, temp = "T", duration = "t",
                       n_total = "n", n_surv = "y",
                       duration_unit = "minutes")
wf <- fit_4pl(wf, chains = 2, iter = 2000, cores = 2, seed = saved$meta$seed,
              refresh = 0, silent = 2, backend = "cmdstanr",
              control = list(adapt_delta = 0.95, max_treedepth = 14))
print(summary(wf$fit))

# Posterior summaries of the four 4PL sub-parameters at T_bar (each is a
# centred-temperature linear model: intercept + slope * (T - T_bar)).
cat("\nPosterior medians of the 4PL sub-parameters at T_bar (centred T):\n")
draws <- posterior::as_draws_df(wf$fit)
post  <- tibble::tibble(
  parameter   = c("low", "up", "k", "mid"),
  intercept   = c(median(0.001 + plogis(draws$b_lowraw_Intercept) * 0.498),
                  median(0.501 + plogis(draws$b_upraw_Intercept) * 0.498),
                  median(exp(draws$b_logk_Intercept)),
                  median(draws$b_mid_Intercept)),
  slope_per_C = c(median(draws$b_lowraw_temp_c),
                  median(draws$b_upraw_temp_c),
                  median(draws$b_logk_temp_c),
                  median(draws$b_mid_temp_c))
)
print(post, digits = 4)

# Extract CTmax / z under both threshold definitions.
et_rel <- suppressMessages(extract_tdt(
  wf, target_surv = "relative", t_ref = 60,
  time_multiplier = 1, ndraws = 1000, lethal = TRUE))
et_abs <- suppressMessages(extract_tdt(
  wf, target_surv = "absolute", t_ref = 60,
  time_multiplier = 1, ndraws = 1000, lethal = TRUE))

joint_z_rel  <- et_rel$z$summary$z_median
joint_c_rel  <- et_rel$CTmax$summary$temp_median
joint_z_abs  <- et_abs$z$summary$z_median
joint_c_abs  <- et_abs$CTmax$summary$temp_median

cat(sprintf("\nJoint 4PL extracted (posterior medians):\n"))
cat(sprintf("  RELATIVE:  z = %.4f   CTmax_1hr = %.4f °C\n",
            joint_z_rel, joint_c_rel))
cat(sprintf("  ABSOLUTE:  z = %.4f   CTmax_1hr = %.4f °C\n",
            joint_z_abs, joint_c_abs))

# ============================================================================
# Step 4 — Two-stage (binomial GLM Stage 1).
# ============================================================================
cat("\n============================================================\n")
cat("Step 4 — Two-stage with binomial GLM Stage 1 (per temperature)\n")
cat("============================================================\n")
log10_lt50_bin <- vapply(sort(unique(data$T)), function(Ti) {
  d_i <- subset(data, T == Ti)
  fit <- glm(cbind(y, n - y) ~ log10(t), family = binomial("logit"),
             data = d_i)
  cat(sprintf("  T = %.1f:  intercept = %+.3f   slope = %+.3f   "
              , Ti, coef(fit)[1], coef(fit)[2]))
  cat(sprintf("log10(LT50) = %.3f\n", -coef(fit)[1] / coef(fit)[2]))
  -coef(fit)[1] / coef(fit)[2]
}, numeric(1))

stage2_bin <- lm(log10_lt50_bin ~ sort(unique(data$T)))
slope_bin  <- coef(stage2_bin)[2]
int_bin    <- coef(stage2_bin)[1]
ts_bin_z   <- -1 / slope_bin
ts_bin_c   <- (log10(60) - int_bin) / slope_bin
cat(sprintf("\nStage 2 OLS:  slope = %+.4f   intercept = %+.4f\n",
            slope_bin, int_bin))
cat(sprintf("Two-stage (binomial) extracted:  z = %.4f   CTmax_1hr = %.4f °C\n",
            ts_bin_z, ts_bin_c))

# ============================================================================
# Step 5 — Two-stage (beta-binomial GLMM Stage 1).
# ============================================================================
cat("\n============================================================\n")
cat("Step 5 — Two-stage with beta-binomial GLMM Stage 1 (per temperature)\n")
cat("============================================================\n")
log10_lt50_bb <- vapply(sort(unique(data$T)), function(Ti) {
  d_i <- subset(data, T == Ti)
  fit <- glmmTMB(cbind(y, n - y) ~ log10(t),
                 family = betabinomial("logit"), data = d_i)
  cat(sprintf("  T = %.1f:  intercept = %+.3f   slope = %+.3f   phi = %.2f   ",
              Ti, fixef(fit)$cond[1], fixef(fit)$cond[2], sigma(fit)))
  cat(sprintf("log10(LT50) = %.3f\n",
              -fixef(fit)$cond[1] / fixef(fit)$cond[2]))
  -fixef(fit)$cond[1] / fixef(fit)$cond[2]
}, numeric(1))

stage2_bb <- lm(log10_lt50_bb ~ sort(unique(data$T)))
slope_bb  <- coef(stage2_bb)[2]
int_bb    <- coef(stage2_bb)[1]
ts_bb_z   <- -1 / slope_bb
ts_bb_c   <- (log10(60) - int_bb) / slope_bb
cat(sprintf("\nStage 2 OLS:  slope = %+.4f   intercept = %+.4f\n",
            slope_bb, int_bb))
cat(sprintf("Two-stage (beta-binomial) extracted:  z = %.4f   CTmax_1hr = %.4f °C\n",
            ts_bb_z, ts_bb_c))

# ============================================================================
# Step 6 — Reproducibility check: regenerated estimates vs the saved row.
#          The simulation default is target_surv = "relative", so the joint
#          row should match the RELATIVE extract.
# ============================================================================
cat("\n============================================================\n")
cat("Step 6 — Reproducibility check (regenerated vs saved row)\n")
cat("============================================================\n")

saved_wide <- saved$row %>%
  dplyr::select(method, quantity, estimate) %>%
  tidyr::pivot_wider(names_from = quantity, values_from = estimate)

reconstructed <- tibble::tribble(
  ~method,                ~z_recon,    ~CTmax_recon,
  "joint_4pl",            joint_z_rel, joint_c_rel,
  "two_stage_bin",        ts_bin_z,    ts_bin_c,
  "two_stage_bb",         ts_bb_z,     ts_bb_c
)

check <- saved_wide %>%
  rename(z_saved = z, CTmax_saved = CTmax_1hr) %>%
  left_join(reconstructed, by = "method") %>%
  mutate(
    dz     = z_recon     - z_saved,
    dCTmax = CTmax_recon - CTmax_saved
  )
print(check, digits = 6)

eps <- 1e-3
ok <- with(check,
           all(abs(dz) < eps, na.rm = TRUE) &&
           all(abs(dCTmax) < eps, na.rm = TRUE))
cat(sprintf("\nReproducibility: %s (tolerance = %g)\n",
            if (ok) "EXACT MATCH" else "MISMATCH — investigate",
            eps))

# ============================================================================
# Step 7 — Side-by-side: every estimator vs analytical truth.
#          For the joint 4PL, both definitions are shown so the
#          relative-vs-absolute gap is visible.
# ============================================================================
cat("\n============================================================\n")
cat("Step 7 — Side-by-side: estimates vs analytical truth\n")
cat("============================================================\n")
out <- tibble::tibble(
  estimator = c("Joint 4PL (relative)", "Joint 4PL (absolute)",
                "Two-stage (binomial)", "Two-stage (beta-binomial)",
                "Truth (absolute)"),
  z         = c(joint_z_rel, joint_z_abs, ts_bin_z, ts_bb_z, truth$z_true),
  CTmax_1hr = c(joint_c_rel, joint_c_abs, ts_bin_c, ts_bb_c,
                truth$CTmax_1hr_true)
) %>%
  mutate(
    z_bias     = z - truth$z_true,
    CTmax_bias = CTmax_1hr - truth$CTmax_1hr_true
  )
print(out, digits = 5)
