# Gated per-group tests for the MIDPOINT parameterisation's `by=` shortcut and
# explicit moderator formulas (the counterpart to test-direct-grouped.R, which
# covers the direct CTmax/z mode). A 2-group fit with a deliberately large
# per-group z gap (A z=5.56, B z=3.33) confirms: (a) `by=` records the moderator
# in meta$group_vars; (b) per-group z/CTmax recover truth via tls(by=); (c) the
# `by=` sugar fits the same model as the four explicit `~ temp_c * grp` formulas.
# Set RUN_BRMS_TESTS=true.

fit_midpoint_by <- function(file_id, ..., seed = 1) {
  raw <- simulate_tdt_2group(seed = 20260622)
  std <- standardize_data(raw, temp = "T", duration = "t_hours",
                          n_total = "n_trials", n_surv = "y_alive",
                          duration_unit = "hours")
  wf <- fit_4pl(std, t_ref = 60, ...,
                chains = 2, iter = 1500, warmup = 750, cores = 2, seed = seed,
                control = list(adapt_delta = 0.95, max_treedepth = 12),
                file = file.path(here::here("tests", "testthat", "fixtures"), file_id),
                refresh = 0)
  attr(wf, "truth_2group") <- attr(raw, "truth_2group")
  wf
}

test_that("midpoint by= records the moderator and recovers per-group truth via tls()", {
  skip_unless_brms()
  wf <- fit_midpoint_by("sim_fit_mid_by_grp", by = "grp")

  expect_equal(wf$meta$parameterization, "midpoint")
  expect_equal(wf$meta$group_vars, "grp")
  expect_true(wf$meta$grouped)
  # all four sub-parameters carry temp_c * grp
  for (p in c("lowraw", "upraw", "logk", "mid"))
    expect_match(deparse(wf$formula$pforms[[p]][[3]]), "temp_c \\* grp")

  tr  <- attr(wf, "truth_2group")
  out <- tls(wf, by = "grp", params = c("z", "ctmax"), t_ref = 60, ndraws = 1000,
             temp_mean = wf$meta$temp_mean)$summary
  gg  <- function(g, q) out$median[out$grp == g & out$quantity == q]
  expect_true("grp" %in% names(out))
  expect_setequal(out$grp, c("A", "B"))

  zA <- gg("A", "z"); zB <- gg("B", "z")
  cA <- gg("A", "CTmax"); cB <- gg("B", "CTmax")
  expect_lt(abs(zA - tr$A$z), 1.5)
  expect_lt(abs(zB - tr$B$z), 1.5)
  expect_gt(zA - zB, 1.0)                                    # large gap -> real mapping
  expect_lt(abs(cA - tr$A$CTmax_1hr), 1.0)
  expect_lt(abs(cB - tr$B$CTmax_1hr), 1.0)
})

test_that("midpoint by= sugar fits the same model as the four explicit formulas", {
  skip_unless_brms()
  wf_by  <- fit_midpoint_by("sim_fit_mid_by_grp", by = "grp")
  wf_exp <- fit_midpoint_by("sim_fit_mid_explicit_grp",
                            low = ~ temp_c * grp, up = ~ temp_c * grp,
                            k   = ~ temp_c * grp, mid = ~ temp_c * grp)
  z_by  <- tls(wf_by,  params = "z", by = "grp", ndraws = NULL,
               temp_mean = wf_by$meta$temp_mean)$summary
  z_exp <- tls(wf_exp, params = "z", by = "grp", ndraws = NULL,
               temp_mean = wf_exp$meta$temp_mean)$summary
  gm <- function(s, g) s$median[s$grp == g & s$quantity == "z"]
  for (gg in c("A", "B"))
    expect_lt(abs(gm(z_by, gg) - gm(z_exp, gg)), 0.3)
})

test_that("explicit low/up/k moderator without mid warns about pooled z", {
  skip_unless_brms()
  raw <- simulate_tdt_2group(seed = 20260622)
  std <- standardize_data(raw, temp = "T", duration = "t_hours",
                          n_total = "n_trials", n_surv = "y_alive",
                          duration_unit = "hours")
  expect_warning(
    fit_4pl(std, low = ~ temp_c * grp, up = ~ temp_c * grp, k = ~ temp_c * grp,
            t_ref = 60, fit = FALSE),
    "mid does not.*POOLED")
})
