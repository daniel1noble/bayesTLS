# Gated per-group tests for the MIDPOINT parameterisation's `by=` shortcut and
# explicit moderator formulas (the counterpart to test-direct-grouped.R, which
# covers the direct CTmax/z mode). A 2-group fit with a deliberately large
# per-group z gap (A z=5.56, B z=3.33) confirms: (a) `by=` records the moderator
# in meta$group_vars so extract_tdt() auto-groups; (b) per-group z/CTmax recover
# truth; (c) the `by=` sugar fits the same model as the four explicit
# `~ temp_c * grp` formulas; (d) extract_tdt(by=) agrees with tls(by=).
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

test_that("midpoint by= records the moderator and auto-groups in extract_tdt()", {
  skip_unless_brms()
  wf <- fit_midpoint_by("sim_fit_mid_by_grp", by = "grp")

  expect_equal(wf$meta$parameterization, "midpoint")
  expect_equal(wf$meta$group_vars, "grp")
  expect_true(wf$meta$grouped)
  # all four sub-parameters carry temp_c * grp
  for (p in c("lowraw", "upraw", "logk", "mid"))
    expect_match(deparse(wf$formula$pforms[[p]][[3]]), "temp_c \\* grp")

  tr  <- attr(wf, "truth_2group")
  out <- extract_tdt(wf, t_ref = 60, ndraws = 1000)          # by = NULL -> auto
  expect_true("grp" %in% names(out$z$summary))
  expect_setequal(out$z$summary$grp, c("A", "B"))

  zA <- out$z$summary$z_median[out$z$summary$grp == "A"]
  zB <- out$z$summary$z_median[out$z$summary$grp == "B"]
  cA <- out$CTmax$summary$temp_median[out$CTmax$summary$grp == "A"]
  cB <- out$CTmax$summary$temp_median[out$CTmax$summary$grp == "B"]
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
  z_by  <- extract_tdt(wf_by,  ndraws = NULL)$z$summary
  z_exp <- extract_tdt(wf_exp, ndraws = NULL)$z$summary
  for (gg in c("A", "B"))
    expect_lt(abs(z_by$z_median[z_by$grp == gg] - z_exp$z_median[z_exp$grp == gg]), 0.3)
})

test_that("midpoint extract_tdt(by=) agrees with tls(by=) per group", {
  skip_unless_brms()
  wf <- fit_midpoint_by("sim_fit_mid_by_grp", by = "grp")
  et <- extract_tdt(wf, t_ref = 60, ndraws = NULL)
  tl <- tls(wf, by = "grp", params = c("z", "ctmax"), mode = "relative",
            t_ref = 60, temp_mean = wf$meta$temp_mean)$summary
  for (gg in c("A", "B")) {
    z_et <- et$z$summary$z_median[et$z$summary$grp == gg]
    z_tl <- tl$median[tl$grp == gg & tl$quantity == "z"]
    expect_lt(abs(z_et - z_tl), 0.3)
  }
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
