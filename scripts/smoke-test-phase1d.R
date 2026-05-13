suppressPackageStartupMessages({
  library(brms); library(dplyr); library(tibble); library(tidyr)
  library(ggplot2); library(utils)
})

lib_files <- c("utils.R", "standardize_data.R", "priors.R", "fit_4pl.R",
               "predict_survival_curves.R", "extract_tdt.R", "tdt_landscape.R",
               "repair.R", "temperature_scenarios.R", "predict_heat_injury.R",
               "plotting.R", "diagnostics.R")
for (f in lib_files) source(file.path("R", f))
cat("sourced:", length(lib_files), "files\n\n")

# Rebuild workflow from cached Phase 1b fit (same as Phase 1c smoke test)
fit <- readRDS(here::here("output", "models", "smoke_phase1b", "fit_smoke1b.rds"))

set.seed(20260512)
truth <- list(ell = 0.02, u = 0.98, k = 6, m_beta0 = 1.5, m_beta1 = -0.18, T_bar = 33)
temps     <- c(30, 32, 34, 36)
durations <- c(0.05, 0.5, 5, 50, 500) / 60
design <- expand.grid(T = temps, t_hours = durations, rep = 1:8) |>
  mutate(log10_t_min = log10(t_hours * 60),
         T_c2        = T - truth$T_bar,
         mid_t       = truth$m_beta0 + truth$m_beta1 * T_c2,
         p_true      = truth$ell + (truth$u - truth$ell) /
                       (1 + exp(truth$k * (log10_t_min - mid_t))),
         n_trials    = 30L,
         p_draw      = rbeta(n(), p_true * 8, (1 - p_true) * 8),
         y_alive     = rbinom(n(), size = n_trials, prob = p_draw))
std <- standardize_data(design, temp = "T", duration = "t_hours",
                        n_total = "n_trials", n_surv = "y_alive")
wf <- structure(list(
  fit = fit, data = std,
  formula = make_4pl_formula(),
  prior   = make_4pl_priors(std),
  meta    = list(
    temp_mean = mean(std$temp), duration_unit = "hours",
    random_effects = NULL, lower = 0, upper = 1,
    bounds = compute_4pl_bounds(0, 1)
  )
), class = "bayes_tls")

fig_dir <- here::here("output", "figs", "smoke_phase1d")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

save_one <- function(p, name) {
  ggsave(file.path(fig_dir, paste0(name, ".png")),
         p, width = 7, height = 4.5, dpi = 120)
}

# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------

cat("=== diagnose_tdt_fit ===\n")
print(diagnose_tdt_fit(wf))

cat("\n=== tdt_parameter_table ===\n")
print(tdt_parameter_table(wf))

# -----------------------------------------------------------------------------
# Plotting: survival curves
# -----------------------------------------------------------------------------

cat("\n=== plot_survival_curves ===\n")
psc <- predict_survival_curves(wf, temps = c(30, 32, 34, 36), ndraws = 300)
p1 <- plot_survival_curves(psc, observed = wf$data)
stopifnot(inherits(p1, "ggplot"))
save_one(p1, "survival_curves")
cat("  ok\n")

# -----------------------------------------------------------------------------
# LTx curve
# -----------------------------------------------------------------------------

cat("=== plot_ltx_curve ===\n")
ltx <- derive_ltx_curve(wf, temp_grid = seq(29, 37, by = 0.5),
                        target_surv = 0.5, ndraws = 300)
p2 <- plot_ltx_curve(ltx)
stopifnot(inherits(p2, "ggplot"))
save_one(p2, "ltx_curve")
cat("  ok\n")

# -----------------------------------------------------------------------------
# TDT landscape
# -----------------------------------------------------------------------------

cat("=== plot_tdt_landscape ===\n")
lsp <- derive_tdt_landscape(wf, ndraws = 200)
p3 <- plot_tdt_landscape(lsp, observed = wf$data)
stopifnot(inherits(p3, "ggplot"))
save_one(p3, "tdt_landscape")
cat("  ok\n")

# -----------------------------------------------------------------------------
# CTmax / T_crit density
# -----------------------------------------------------------------------------

cat("=== plot_temperature_density (CTmax + T_crit) ===\n")
et <- extract_tdt(wf, ndraws = 300)
p4 <- plot_temperature_density(et$CTmax,
                                x_label = expression(CT[max[1*hr]]~("°C")))
p5 <- plot_temperature_density(et$T_crit,
                                x_label = expression(T[crit]~("°C")))
stopifnot(inherits(p4, "ggplot"), inherits(p5, "ggplot"))
save_one(p4, "ctmax_density")
save_one(p5, "tcrit_density")
cat("  ok\n")

# -----------------------------------------------------------------------------
# Temperature scenarios + heat injury trajectories
# -----------------------------------------------------------------------------

cat("=== plot_temperature_scenarios ===\n")
scens <- make_temperature_scenarios(baseline = 20, spike_temp = 28,
                                    n_hours = 96,
                                    spike_times_single = 24,
                                    spike_times_multi = c(24, 48, 72))
p6 <- plot_temperature_scenarios(scens, T_c = 24)
stopifnot(inherits(p6, "ggplot"))
save_one(p6, "scenarios")
cat("  ok\n")

cat("=== plot_heat_injury (single_spike) ===\n")
hi_single <- predict_heat_injury(scens$single_spike, wf, T_c = 24, ndraws = 200)
p7 <- plot_heat_injury(hi_single)
save_one(p7, "heat_injury_single")
cat("  ok\n")

cat("=== plot_heat_injury (multi_spike with repair) ===\n")
repair_pars <- list(TA = 14065, TAL = 50000, TAH = 120000,
                    TL = 10.5 + 273.15, TH = 22.5 + 273.15,
                    TREF = 17 + 273.15, r_ref = 0.005)
hi_multi_rp <- predict_heat_injury(scens$multi_spike, wf, T_c = 24,
                                   ndraws = 200, repair = TRUE,
                                   repair_pars = repair_pars)
p8 <- plot_heat_injury(hi_multi_rp)
save_one(p8, "heat_injury_multi_repair")
cat("  ok\n")

cat("=== plot_repair_rate ===\n")
p9 <- plot_repair_rate(seq(5, 35, length.out = 200), repair_pars)
stopifnot(inherits(p9, "ggplot"))
save_one(p9, "repair_rate")
cat("  ok\n")

cat("\nfigures saved to:", fig_dir, "\n")
cat(length(list.files(fig_dir, pattern = "\\.png$")), "png files\n")
cat("\nALL CHECKS COMPLETE\n")
