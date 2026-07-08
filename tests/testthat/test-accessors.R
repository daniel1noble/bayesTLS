# Tests for the kept accessors. Input-validation tests are brms-free; the
# round-trip tests against predict_*() are gated behind RUN_BRMS_TESTS.
# get_tls_est() reads a `tls` object (the result of tls()); get_4pl_est() reads
# a fitted workflow; get_hi_draws()/get_surv_draws() read predict_*() results.

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

test_that("get_tls_est pulls summary/draws and filters params from a tls object", {
  fake <- structure(
    list(
      summary = tibble::tibble(
        quantity = c("z", "CTmax", "Tcrit"),
        median   = c(3, 40, 30), lower = c(2, 39, 29), upper = c(4, 41, 31)),
      draws = tibble::tibble(
        quantity = rep(c("z", "CTmax", "Tcrit"), each = 5),
        .draw    = rep(1:5, 3), value = as.numeric(1:15)),
      meta = list()),
    class = c("tls", "list"))

  s <- get_tls_est(fake, "summary")
  expect_true(all(c("quantity", "median", "lower", "upper") %in% names(s)))
  expect_equal(nrow(s), 3L)
  expect_s3_class(s, "tbl_df")

  d <- get_tls_est(fake, "draws")
  expect_true(all(c("quantity", ".draw", "value") %in% names(d)))
  expect_equal(nrow(d), 15L)

  # default is summary
  expect_equal(get_tls_est(fake), s)

  # params filter, case-insensitive
  expect_equal(unique(get_tls_est(fake, "summary", "z")$quantity), "z")
  expect_setequal(unique(get_tls_est(fake, "draws", c("Z", "ctmax"))$quantity),
                  c("z", "CTmax"))

  # errors
  expect_error(get_tls_est(list(z = 1, CTmax = 2)), "must be a `tls` object")
  expect_error(get_tls_est(fake, "summary", "nope"), "No matching")
})


test_that("get_tls_est passes grouped columns through and filters per group", {
  g <- structure(list(
    summary = tibble::tibble(species = rep(c("A", "B"), each = 3),
       quantity = rep(c("z", "CTmax", "Tcrit"), 2),
       median = 1:6, lower = 0:5, upper = 2:7),
    draws = tibble::tibble(species = rep(c("A", "B"), each = 4),
       quantity = rep(c("z", "CTmax"), 4), .draw = rep(1:2, 4), value = as.numeric(1:8)),
    meta = list()), class = c("tls", "list"))
  expect_true("species" %in% names(get_tls_est(g, "summary")))
  expect_equal(nrow(get_tls_est(g, "summary")), 6L)
  expect_equal(nrow(get_tls_est(g, "summary", "CTmax")), 2L)   # one CTmax row per group
  expect_true("species" %in% names(get_tls_est(g, "draws", "z")))
})

test_that("get_4pl_est validates `what` and requires a fitted workflow", {
  expect_error(get_4pl_est(list(fit = NULL), "nope"))               # match.arg
  expect_error(get_4pl_est(list(fit = NULL), "summary"), "fit")     # has_fit guard
  expect_error(get_4pl_est(list(fit = NULL), "draws"), "fit")
})
