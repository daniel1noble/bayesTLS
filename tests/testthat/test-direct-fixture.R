# Gated brms integration test for the DIRECT CTmax/z parameterisation. Fits a
# small cached single-group direct model on the same simulate_tdt() data as the
# midpoint fixture, then asserts that (a) the rewired coefficient-read helpers
# recover the planted z/CTmax, (b) they agree with the midpoint fixture (the
# equivalence claim of the redesign), and (c) extract_tdt/tdt_parameter_table/
# predict_heat_injury return finite quantities on a direct fit (the integration-
# level version of the silent-NA regression guard). Set RUN_BRMS_TESTS=true.

test_that("extract_tdt recovers z and CTmax_1hr from truth on a DIRECT fit", {
  skip_unless_brms()

  wf  <- load_fixture_workflow_direct()
  ts  <- truth_summary()
  out <- extract_tdt(wf, t_ref = 60, ndraws = 1500)

  expect_lt(abs(out$z$summary$z_median        - ts$z),         1.5)
  expect_lt(abs(out$CTmax$summary$temp_median - ts$CTmax_1hr), 1.0)
  expect_lte(out$z$summary$z_lower,        ts$z)
  expect_gte(out$z$summary$z_upper,        ts$z)
  expect_lte(out$CTmax$summary$temp_lower, ts$CTmax_1hr)
  expect_gte(out$CTmax$summary$temp_upper, ts$CTmax_1hr)
})

test_that("a direct fit returns FINITE z/CTmax/T_crit (integration silent-NA guard)", {
  skip_unless_brms()

  wf  <- load_fixture_workflow_direct()
  out <- suppressMessages(extract_tdt(wf, t_ref = 60, ndraws = 1500, lethal = TRUE))
  # Before the rewire, a direct fit had no b_mid_* -> every quantity was NA.
  expect_true(is.finite(out$z$summary$z_median))
  expect_true(is.finite(out$CTmax$summary$temp_median))
  expect_true(is.finite(out$T_crit$summary$temp_median))
  expect_true(is.finite(out$z$summary$z_lower) && is.finite(out$z$summary$z_upper))
  expect_lt(out$T_crit$summary$temp_median, out$CTmax$summary$temp_median)
})

test_that("direct and midpoint parameterisations agree (equivalence claim)", {
  skip_unless_brms()

  wf_d <- load_fixture_workflow_direct()
  wf_m <- load_fixture_workflow()
  ed <- extract_tdt(wf_d, t_ref = 60, ndraws = 1500)
  em <- extract_tdt(wf_m, t_ref = 60, ndraws = 1500)

  # Same data, two parameterisations -> z and CTmax medians should coincide up
  # to MCMC noise on this small fit. Tolerances are far tighter than the truth-
  # recovery window, so they detect a genuine parameterisation discrepancy.
  expect_lt(abs(ed$z$summary$z_median        - em$z$summary$z_median),        0.6)
  expect_lt(abs(ed$CTmax$summary$temp_median - em$CTmax$summary$temp_median), 0.5)
})

test_that("tdt_parameter_table on a direct fit is consistent with extract_tdt", {
  skip_unless_brms()

  wf  <- load_fixture_workflow_direct()
  tab <- tdt_parameter_table(wf)
  out <- extract_tdt(wf, t_ref = 60, ndraws = 1500)

  expect_true(all(is.finite(tab$median)))
  expect_setequal(
    tab$parameter,
    c("low (lower asymptote)", "up (upper asymptote)", "k (slope)",
      "CTmax (°C, at reference dose)", "z (°C)")
  )
  # z from the table (= exp(logz), at temp_mean) and from extract_tdt agree.
  z_tab <- tab$median[tab$parameter == "z (°C)"]
  expect_lt(abs(z_tab - out$z$summary$z_median), 0.6)
})

test_that("seed makes extract_tdt and tls reproducible on a direct fit (D-T2)", {
  skip_unless_brms()
  wf <- load_fixture_workflow_direct()

  # extract_tdt: ndraws < fixture draws -> subsample + T_crit runif both seeded.
  a  <- suppressMessages(extract_tdt(wf, ndraws = 300, lethal = TRUE, seed = 5))
  b  <- suppressMessages(extract_tdt(wf, ndraws = 300, lethal = TRUE, seed = 5))
  cc <- suppressMessages(extract_tdt(wf, ndraws = 300, lethal = TRUE, seed = 6))
  expect_identical(a$z$summary$z_median,          b$z$summary$z_median)
  expect_identical(a$T_crit$summary$temp_median,  b$T_crit$summary$temp_median)
  expect_false(isTRUE(all.equal(a$T_crit$summary$temp_median,
                                cc$T_crit$summary$temp_median)))

  # tls(): posterior_linpred subsample + T_crit runif both seeded.
  t1 <- tls(wf, lethal = TRUE, ndraws = 300, seed = 5)$summary
  t2 <- tls(wf, lethal = TRUE, ndraws = 300, seed = 5)$summary
  expect_identical(t1$median, t2$median)
})

test_that("derive_z and derive_temperature_for_duration are wired for direct fits", {
  skip_unless_brms()

  wf  <- load_fixture_workflow_direct()
  ts  <- truth_summary()
  out <- extract_tdt(wf, t_ref = 60, ndraws = 1500)

  # derive_z(): direct fit -> finite z recovering exp(logz); agrees with extract_tdt.
  zr <- derive_z(wf)
  expect_true(is.finite(zr$summary$z_median))
  expect_lt(abs(zr$summary$z_median - ts$z), 1.5)
  expect_lt(abs(zr$summary$z_median - out$z$summary$z_median), 0.3)

  # derive_temperature_for_duration(): finite CTmax at the 1 h reference.
  tg <- seq(min(wf$data$temp) - 2, max(wf$data$temp) + 2, by = 0.05)
  dd <- derive_temperature_for_duration(wf, exposure_duration = 1, temp_grid = tg,
                                        target_surv = "relative")
  expect_true(is.finite(dd$summary$temp_median))
  expect_lt(abs(dd$summary$temp_median - out$CTmax$summary$temp_median), 1.0)
})

test_that("predict_heat_injury runs on a direct fit and agrees with midpoint", {
  skip_unless_brms()

  wf_d  <- load_fixture_workflow_direct()
  wf_m  <- load_fixture_workflow()
  trace <- data.frame(time = 0:8,
                      temp = c(28, 30, 32, 34, 35, 34, 32, 30, 28))

  hd <- suppressMessages(predict_heat_injury(trace, wf_d, ndraws = 1500))$summary
  hm <- suppressMessages(predict_heat_injury(trace, wf_m, ndraws = 1500))$summary

  expect_equal(nrow(hd), nrow(trace))          # non-empty (was 0-row before fix)
  expect_true(all(is.finite(hd$surv_median)))
  # final predicted survival agrees between parameterisations
  expect_lt(abs(tail(hd$surv_median, 1) - tail(hm$surv_median, 1)), 0.1)
})
