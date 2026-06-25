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
  expect_error(get_tls_draws(list(foo = 1)),     "extract_tdt")
  expect_error(get_tls_summary(list(foo = 1)),   "extract_tdt")
})

test_that("get_tls_summary reshapes per-quantity summaries to a tidy long table", {
  # Harmonises the quantity-specific column names (z_median vs temp_median) into
  # quantity / median / lower / upper, drops the T_crit rate-floor columns, and
  # labels quantities z / CTmax / Tcrit (matching tls()).
  fake <- list(
    z      = list(summary = tibble::tibble(z_median = 5, z_lower = 4, z_upper = 6)),
    CTmax  = list(summary = tibble::tibble(temp_lower = 30, temp_median = 31, temp_upper = 32)),
    T_crit = list(summary = tibble::tibble(TC_rate_low = 0.1, TC_rate_high = 1,
                                           temp_lower = 20, temp_median = 22, temp_upper = 24))
  )
  s <- get_tls_summary(fake)
  expect_named(s, c("quantity", "median", "lower", "upper"))
  expect_equal(s$quantity, c("z", "CTmax", "Tcrit"))
  expect_equal(s$median,   c(5, 31, 22))
  expect_equal(s$lower,    c(4, 30, 20))
  expect_equal(s$upper,    c(6, 32, 24))
  expect_false(any(grepl("TC_rate", names(s))))      # rate-floor columns dropped

  fake$T_crit <- NULL                                # lethal = FALSE
  expect_equal(get_tls_summary(fake)$quantity, c("z", "CTmax"))
})

test_that("get_tls_draws inner-joins on .draw: keeps only shared draws, never mis-pairs", {
  # Quantities are filtered to finite values independently inside extract_tdt(),
  # so their .draw sets can differ. Here z has draws 1:4, CTmax 1:3 (draw 4
  # dropped), T_crit 2:3. get_tls_draws() must return only the draws present in
  # ALL THREE, with each value carried from its own draw (a column-bind would
  # mis-pair, e.g. give z = 5 against CTmax = 30 on different draws).
  fake <- list(
    z      = list(draws = tibble::tibble(.draw = 1:4, z    = c(5, 6, 7, 8))),
    CTmax  = list(draws = tibble::tibble(.draw = 1:3, temp = c(30, 31, 32))),
    T_crit = list(draws = tibble::tibble(.draw = 2:3, temp = c(20, 21),
                                         log10_rate = c(-2.5, -2.5)))
  )
  d <- get_tls_draws(fake)
  expect_named(d, c(".draw", "z", "CTmax", "T_crit"))
  expect_equal(d$.draw,  2:3)          # intersection of all three draw sets
  expect_equal(d$z,      c(6, 7))      # z FROM draws 2,3 (not 5,6 -> would be mis-paired)
  expect_equal(d$CTmax,  c(31, 32))
  expect_equal(d$T_crit, c(20, 21))
  expect_false("log10_rate" %in% names(d))   # auxiliary column dropped
})

test_that("get_tls_draws omits T_crit when extract_tdt was lethal = FALSE", {
  fake <- list(
    z      = list(draws = tibble::tibble(.draw = 1:3, z    = c(5, 6, 7))),
    CTmax  = list(draws = tibble::tibble(.draw = 1:3, temp = c(30, 31, 32))),
    T_crit = NULL
  )
  d <- get_tls_draws(fake)
  expect_named(d, c(".draw", "z", "CTmax"))
  expect_false("T_crit" %in% names(d))
  expect_equal(nrow(d), 3L)
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

test_that("get_tls_draws merges z/CTmax/T_crit on the SAME draw (joint pairing preserved)", {
  skip_unless_brms()

  wf <- load_fixture_workflow()
  et <- suppressMessages(extract_tdt(wf, t_ref = 60, ndraws = 300, lethal = TRUE))
  d  <- get_tls_draws(et)

  expect_s3_class(d, "tbl_df")
  expect_named(d, c(".draw", "z", "CTmax", "T_crit"))
  expect_false(any(duplicated(d$.draw)))               # one row per draw

  # Pairing: each column's value matches the per-quantity accessor for that draw.
  zd <- get_z_draws(et); cd <- get_ctmax_draws(et); td <- get_tcrit_draws(et)
  expect_equal(d$z,      zd$z[match(d$.draw, zd$.draw)])
  expect_equal(d$CTmax,  cd$CTmax[match(d$.draw, cd$.draw)])
  expect_equal(d$T_crit, td$T_crit[match(d$.draw, td$.draw)])

  # lethal = FALSE -> z and CTmax only, no T_crit column.
  e2 <- extract_tdt(wf, t_ref = 60, ndraws = 300, lethal = FALSE)
  expect_named(get_tls_draws(e2), c(".draw", "z", "CTmax"))
})

test_that("get_tls_summary matches the individual summaries and carries the moderator", {
  skip_unless_brms()

  wf <- load_fixture_workflow()
  et <- suppressMessages(extract_tdt(wf, t_ref = 60, ndraws = 300, lethal = TRUE))
  s  <- get_tls_summary(et)
  expect_named(s, c("quantity", "median", "lower", "upper"))
  expect_equal(s$quantity, c("z", "CTmax", "Tcrit"))
  expect_equal(s$median[s$quantity == "z"],     get_z_summary(et)$z_median)
  expect_equal(s$median[s$quantity == "CTmax"], get_ctmax_summary(et)$temp_median)
  expect_equal(s$upper[s$quantity == "Tcrit"],  get_tcrit_summary(et)$temp_upper)

  # lethal = FALSE -> z and CTmax only
  e2 <- extract_tdt(wf, t_ref = 60, ndraws = 300)
  expect_equal(get_tls_summary(e2)$quantity, c("z", "CTmax"))

  # grouped: one row per quantity x group, moderator preserved
  wg <- load_fixture_workflow_grouped()
  sg <- get_tls_summary(suppressMessages(
    extract_tdt(wg, t_ref = 60, lethal = TRUE, ndraws = NULL)))
  expect_true(all(c("grp", "quantity", "median", "lower", "upper") %in% names(sg)))
  expect_setequal(unique(sg$grp), c("A", "B"))
  expect_equal(nrow(sg), 6L)                          # 3 quantities x 2 groups
})

test_that("get_tls_draws keeps the moderator and joins WITHIN group (no cross-join)", {
  skip_unless_brms()

  wf <- load_fixture_workflow_grouped()
  et <- suppressMessages(extract_tdt(wf, t_ref = 60, lethal = TRUE, ndraws = NULL))
  d  <- get_tls_draws(et)

  expect_true(all(c("grp", ".draw", "z", "CTmax", "T_crit") %in% names(d)))
  expect_setequal(unique(d$grp), c("A", "B"))
  # No duplicate (grp, .draw) keys: a many-to-many join would scramble groups
  # and inflate rows to ~n^2; the correct inner join is 1:1 per (grp, .draw).
  expect_false(any(duplicated(d[, c("grp", ".draw")])))
  expect_lte(nrow(d), nrow(get_z_draws(et)))
  # within group A, z stays paired with its own draw
  da <- d[d$grp == "A", ]; za <- get_z_draws(et); za <- za[za$grp == "A", ]
  expect_equal(da$z[order(da$.draw)], za$z[order(za$.draw)])
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

test_that("grouped extract_tdt accessors preserve the moderator column (key on grp,.draw)", {
  skip_unless_brms()
  wf <- load_fixture_workflow_grouped()
  et <- extract_tdt(wf, t_ref = 60, lethal = TRUE, ndraws = NULL)
  for (acc in list(get_z_draws, get_ctmax_draws, get_tcrit_draws)) {
    d <- acc(et)
    expect_true("grp" %in% names(d))                 # group column carried, not dropped
    expect_setequal(unique(d$grp), c("A", "B"))
    expect_gt(sum(duplicated(d$.draw)), 0L)           # .draw repeats across groups
  }
  # The headline use case: a per-group z-vs-CTmax contrast keyed on (grp, .draw)
  # must be 1:1 (the bug produced a many-to-many join that scrambled groups).
  zj <- merge(get_z_draws(et), get_ctmax_draws(et), by = c("grp", ".draw"))
  expect_equal(nrow(zj), nrow(get_z_draws(et)))
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
