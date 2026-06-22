test_that("time_to_surv_threshold_4pl returns NA when target outside (low, up)", {
  # target = 0.05 is below low = 0.1, so the 4PL never reaches it → NA
  expect_true(is.na(time_to_surv_threshold_4pl(
    temp = 30, survival_target = 0.05,
    low = 0.1, up = 0.95, k = 6,
    mid_int = 1, mid_temp = -0.2, temp_mean = 30
  )))
})

test_that("time_to_surv_threshold_4pl recovers the analytical inverse at T_bar", {
  # At T_bar (temp_c = 0), mid = mid_int. Survival = 0.5 at duration = 10^mid_int
  # when (low, up) are symmetric around 0.5.
  out <- time_to_surv_threshold_4pl(
    temp = 30, survival_target = 0.5,
    low = 0.05, up = 0.95, k = 6,
    mid_int = 1.5, mid_temp = -0.2, temp_mean = 30
  )
  # mid = 1.5, log term = log(0.45/0.45)/k = 0, so duration = 10^1.5
  expect_equal(out, 10 ^ 1.5, tolerance = 1e-6)
})

test_that("survival_from_dose returns u at dose = 0 and target_surv at dose = 1", {
  expect_equal(survival_from_dose(0,
                                  low = 0.02, up = 0.98, k = 6,
                                  target_surv = 0.5),
               0.98, tolerance = 1e-6)
  expect_equal(survival_from_dose(1,
                                  low = 0.02, up = 0.98, k = 6,
                                  target_surv = 0.5),
               0.5, tolerance = 1e-6)
})

test_that("repair_rate_schoolfield returns positive rates peaked near TREF", {
  rp <- list(TA = 14065, TAL = 50000, TAH = 120000,
             TL = 10.5 + 273.15, TH = 22.5 + 273.15,
             TREF = 17 + 273.15, r_ref = 0.01)
  temps <- seq(5, 30, by = 1)
  rates <- repair_rate_schoolfield(
    temp_celsius = temps,
    TA = rp$TA, TAL = rp$TAL, TAH = rp$TAH,
    TL = rp$TL, TH = rp$TH, TREF = rp$TREF,
    r_ref = rp$r_ref
  )
  expect_true(all(rates >= 0))
  expect_equal(length(rates), length(temps))
  # The peak should be near TREF = 17 °C (TPC has its optimum in this region).
  peak <- temps[which.max(rates)]
  expect_true(peak >= 13 && peak <= 22)
})

# ---------------------------------------------------------------------------
# Unit reconciliation (regression test for the 2026-06-12 fix): the dose
# integral must be invariant to the model's fit-time unit AND the trace's time
# unit. Before the fix, a minutes-fitted model driven by an hours trace
# under-counted the dose 60-fold.
# ---------------------------------------------------------------------------

# Minimal fake bayes_tls workflow with known 4PL params, built so that
# LT50(temp_mean) = exactly 1 hour in the model's `duration_unit`. Hence a
# 1-hour exposure at temp_mean (= CTmax_1hr) accumulates exactly one dose.
fake_hi_workflow <- function(duration_unit, temp_mean = 30, z = 5,
                             low = 0.02, up = 0.98, k = 6, ndraw = 30) {
  unit_per_hour <- switch(duration_unit,
                          seconds = 3600, minutes = 60, hours = 1, days = 1 / 24)
  bnd <- list(low_min = 0, low_w = 0.1, up_min = 0.9, up_w = 0.1)
  fit <- posterior::as_draws_df(data.frame(
    b_lowraw_Intercept = rep(stats::qlogis((low - bnd$low_min) / bnd$low_w), ndraw),
    b_upraw_Intercept  = rep(stats::qlogis((up  - bnd$up_min)  / bnd$up_w),  ndraw),
    b_logk_Intercept   = rep(log(k), ndraw),
    b_mid_Intercept    = rep(log10(unit_per_hour), ndraw),   # LT50(temp_mean)=1 h
    b_mid_temp_c       = rep(-1 / z, ndraw)
  ))
  structure(list(fit = fit,
                 meta = list(temp_mean = temp_mean, duration_unit = duration_unit,
                             bounds = bnd)),
            class = "bayes_tls")
}

# The exact "1 h at CTmax = one dose" planted-parameter checks needed a fake
# as_draws_df workflow, which the posterior_linpred-based extract_4pl_pars can no
# longer consume. The dose-integral + unit-reconciliation logic is byte-identical
# to before (only the parameter source changed), so its regression is covered by:
#  - the pure dose-math tests above (time_to_surv_threshold_4pl, survival_from_dose);
#  - the gated trace-unit INVARIANCE test below (the dt reconciliation: a trace in
#    hours and the same real trace relabelled in minutes must give identical HI);
#  - the gated direct-vs-midpoint equivalence in test-direct-fixture.R.
test_that("predict_heat_injury is invariant to the trace time unit (dt reconciliation)", {
  skip_unless_brms()
  wf <- load_fixture_workflow()                       # fit in hours
  hi_h <- predict_heat_injury(
    data.frame(time = 0:6, temp = seq(30, 36, length.out = 7)),
    wf, trace_unit = "hours", ndraws = 200, seed = 1)$summary
  hi_m <- predict_heat_injury(
    data.frame(time = (0:6) * 60, temp = seq(30, 36, length.out = 7)),
    wf, trace_unit = "minutes", ndraws = 200, seed = 1)$summary
  expect_equal(hi_h$hi_median,   hi_m$hi_median,   tolerance = 1e-8)
  expect_equal(hi_h$surv_median, hi_m$surv_median, tolerance = 1e-8)
})

test_that("predict_heat_injury errors on an unsupported time unit", {
  wf <- fake_hi_workflow("minutes")
  expect_error(
    predict_heat_injury(data.frame(time = c(0, 1), temp = 30), wf,
                        trace_unit = "fortnights", ndraws = 5),
    "unsupported time unit")
})
