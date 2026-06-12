# Tests for the get_*_draws() accessors. Input-validation tests are
# brms-free; the round-trip tests against extract_tdt() / predict_*() are
# gated behind RUN_BRMS_TESTS.

test_that("get_*_draws / get_*_summary reject non-extract_tdt input", {
  expect_error(get_z_draws(list(foo = 1)),       "extract_tdt")
  expect_error(get_ctmax_draws(list(foo = 1)),   "extract_tdt")
  expect_error(get_tcrit_draws(list(foo = 1)),   "extract_tdt")
  expect_error(get_z_summary(list(foo = 1)),     "extract_tdt")
  expect_error(get_ctmax_summary(list(foo = 1)), "extract_tdt")
  expect_error(get_tcrit_summary(list(foo = 1)), "extract_tdt")
})

test_that("get_tcrit_draws / get_tcrit_summary error when T_crit is absent (lethal = FALSE)", {
  fake <- list(z     = list(draws = tibble::tibble(.draw = 1, z = 1.0)),
               CTmax = list(draws = tibble::tibble(.draw = 1, temp = 30)),
               T_crit = NULL)
  expect_error(get_tcrit_draws(fake),   "lethal = TRUE")
  expect_error(get_tcrit_summary(fake), "lethal = TRUE")
})

test_that("get_hi_draws errors helpfully when save_draws was FALSE", {
  fake <- list(
    summary = tibble::tibble(time = c(0, 1), temp = c(20, 20),
                             hi_median = c(0, 0)),
    meta    = list(),
    draws   = NULL
  )
  expect_error(get_hi_draws(fake), "save_draws = TRUE")
})

test_that("get_surv_draws dispatches on input shape", {
  expect_error(get_surv_draws(list(foo = 1)),
               "predict_survival_curves|predict_heat_injury")
})

test_that("get_z_draws and get_ctmax_draws round-trip with extract_tdt", {
  skip_unless_brms()

  wf <- load_fixture_workflow()
  et <- extract_tdt(wf, t_ref = 60, ndraws = 300)

  zd <- get_z_draws(et)
  expect_s3_class(zd, "tbl_df")
  expect_named(zd, c(".draw", "z"))
  expect_equal(nrow(zd), nrow(et$z$draws))

  cd <- get_ctmax_draws(et)
  expect_named(cd, c(".draw", "CTmax"))
  expect_equal(nrow(cd), nrow(et$CTmax$draws))

  zs <- get_z_summary(et)
  expect_s3_class(zs, "tbl_df")
  expect_named(zs, c("z_median", "z_lower", "z_upper"))
  expect_equal(nrow(zs), 1L)

  cs <- get_ctmax_summary(et)
  expect_true(all(c("temp_lower", "temp_median", "temp_upper") %in% names(cs)))
  expect_gte(nrow(cs), 1L)
})

test_that("get_tcrit_draws round-trips with extract_tdt(lethal = TRUE)", {
  skip_unless_brms()

  wf <- load_fixture_workflow()
  et <- suppressMessages(extract_tdt(wf, t_ref = 60, ndraws = 300,
                                     lethal = TRUE))
  td <- get_tcrit_draws(et)
  expect_named(td, c(".draw", "T_crit", "log10_rate"))
  expect_equal(nrow(td), nrow(et$T_crit$draws))

  ts <- get_tcrit_summary(et)
  expect_named(ts, c("TC_rate_low", "TC_rate_high",
                     "temp_lower", "temp_median", "temp_upper"))
  expect_equal(nrow(ts), 1L)
})

test_that("get_hi_draws and get_surv_draws round-trip with predict_heat_injury", {
  skip_unless_brms()

  wf    <- load_fixture_workflow()
  scens <- make_temperature_scenarios(baseline = 20, spike_temp = 28,
                                      n_hours = 24,
                                      spike_times_single = 12,
                                      spike_times_multi  = c(12, 18))
  hi <- predict_heat_injury(scens$single_spike, wf,
                            T_c = 24, ndraws = 100, save_draws = TRUE)

  hd <- get_hi_draws(hi)
  expect_true(all(c(".draw", "time", "temp", "hi", "survival") %in%
                  names(hd)))
  expect_gt(nrow(hd), 0)

  sd <- get_surv_draws(hi)
  expect_named(sd, c(".draw", "time", "temp", "survival"))
  expect_equal(nrow(sd), nrow(hd))
})

test_that("get_surv_draws round-trips with predict_survival_curves", {
  skip_unless_brms()

  wf  <- load_fixture_workflow()
  psc <- predict_survival_curves(wf,
                                 temps     = c(32, 34),
                                 durations = c(0.5, 1),
                                 ndraws    = 100)

  sd <- get_surv_draws(psc)
  expect_named(sd, c(".draw", "temp", "duration", "survival"))
  expect_equal(nrow(sd), nrow(psc$grid) * 100)
})
