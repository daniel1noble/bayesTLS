# Gated brms integration test for the DIRECT CTmax/z parameterisation. Fits a
# small cached single-group direct model on the same simulate_tdt() data as the
# midpoint fixture, then asserts that (a) the rewired coefficient-read helpers
# recover the planted z/CTmax, (b) they agree with the midpoint fixture (the
# equivalence claim of the redesign), and (c) tls/tdt_parameter_table/
# predict_heat_injury return finite quantities on a direct fit (the integration-
# level version of the silent-NA regression guard). Set RUN_BRMS_TESTS=true.

test_that("tls() recovers z and CTmax_1hr from truth on a DIRECT fit", {
  skip_unless_brms()

  wf  <- load_fixture_workflow_direct()
  ts  <- truth_summary()
  out <- tls(wf, params = c("z", "ctmax"), t_ref = 60, ndraws = 1500)$summary
  g   <- function(q, col) out[[col]][out$quantity == q]

  expect_lt(abs(g("z", "median")     - ts$z),         1.5)
  expect_lt(abs(g("CTmax", "median") - ts$CTmax_1hr), 1.0)
  expect_lte(g("z", "lower"),     ts$z)
  expect_gte(g("z", "upper"),     ts$z)
  expect_lte(g("CTmax", "lower"), ts$CTmax_1hr)
  expect_gte(g("CTmax", "upper"), ts$CTmax_1hr)
})

test_that("a direct fit returns FINITE z/CTmax/T_crit (integration silent-NA guard)", {
  skip_unless_brms()

  wf  <- load_fixture_workflow_direct()
  out <- suppressMessages(
    tls(wf, params = "all", t_ref = 60, ndraws = 1500, lethal = TRUE))$summary
  g   <- function(q, col) out[[col]][out$quantity == q]
  # Before the rewire, a direct fit had no b_mid_* -> every quantity was NA.
  expect_true(is.finite(g("z", "median")))
  expect_true(is.finite(g("CTmax", "median")))
  expect_true(is.finite(g("Tcrit", "median")))
  expect_true(is.finite(g("z", "lower")) && is.finite(g("z", "upper")))
  expect_lt(g("Tcrit", "median"), g("CTmax", "median"))
})

test_that("direct and midpoint parameterisations agree (equivalence claim)", {
  skip_unless_brms()

  wf_d <- load_fixture_workflow_direct()
  wf_m <- load_fixture_workflow()
  ed <- tls(wf_d, params = c("z", "ctmax"), t_ref = 60, ndraws = 1500)$summary
  em <- tls(wf_m, params = c("z", "ctmax"), t_ref = 60, ndraws = 1500)$summary
  g  <- function(s, q) s$median[s$quantity == q]

  # Same data, two parameterisations -> z and CTmax medians should coincide up
  # to MCMC noise on this small fit. Tolerances are far tighter than the truth-
  # recovery window, so they detect a genuine parameterisation discrepancy.
  expect_lt(abs(g(ed, "z")     - g(em, "z")),     0.6)
  expect_lt(abs(g(ed, "CTmax") - g(em, "CTmax")), 0.5)
})

test_that("tdt_parameter_table on a direct fit is consistent with tls", {
  skip_unless_brms()

  wf  <- load_fixture_workflow_direct()
  tab <- tdt_parameter_table(wf)
  out <- tls(wf, params = "z", t_ref = 60, ndraws = 1500)$summary

  expect_true(all(is.finite(tab$median)))
  expect_setequal(
    tab$parameter,
    c("low (lower asymptote)", "up (upper asymptote)", "k (slope)",
      "CTmax (°C, at reference dose)", "z (°C)")
  )
  # z from the table (= exp(logz), at temp_mean) and from tls() agree.
  z_tab <- tab$median[tab$parameter == "z (°C)"]
  expect_lt(abs(z_tab - out$median[out$quantity == "z"]), 0.6)
})

test_that("seed makes tls() reproducible on a direct fit (D-T2)", {
  skip_unless_brms()
  wf <- load_fixture_workflow_direct()

  # tls(): posterior_linpred subsample + T_crit runif both seeded, so a fixed
  # seed reproduces exactly and a different seed perturbs the derived quantities.
  t1 <- tls(wf, lethal = TRUE, ndraws = 300, seed = 5)$summary
  t2 <- tls(wf, lethal = TRUE, ndraws = 300, seed = 5)$summary
  t3 <- tls(wf, lethal = TRUE, ndraws = 300, seed = 6)$summary
  gt <- function(s, q) s$median[s$quantity == q]
  expect_identical(t1$median, t2$median)
  expect_identical(gt(t1, "z"), gt(t2, "z"))
  expect_false(isTRUE(all.equal(gt(t1, "Tcrit"), gt(t3, "Tcrit"))))
})

test_that("derive_z and derive_temperature_for_duration are wired for direct fits", {
  skip_unless_brms()

  wf  <- load_fixture_workflow_direct()
  ts  <- truth_summary()
  out <- tls(wf, params = c("z", "ctmax"), t_ref = 60, ndraws = 1500)$summary
  g   <- function(q) out$median[out$quantity == q]

  # derive_z(): direct fit -> finite z recovering exp(logz); agrees with tls().
  zr <- derive_z(wf)
  expect_true(is.finite(zr$summary$z_median))
  expect_lt(abs(zr$summary$z_median - ts$z), 1.5)
  expect_lt(abs(zr$summary$z_median - g("z")), 0.3)

  # derive_temperature_for_duration(): finite CTmax at the 1 h reference.
  tg <- seq(min(wf$data$temp) - 2, max(wf$data$temp) + 2, by = 0.05)
  dd <- derive_temperature_for_duration(wf, exposure_duration = 1, temp_grid = tg,
                                        target_surv = "relative")
  expect_true(is.finite(dd$summary$temp_median))
  expect_lt(abs(dd$summary$temp_median - g("CTmax")), 1.0)
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
