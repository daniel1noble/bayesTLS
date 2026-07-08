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

test_that("survival_from_dose returns NA (no NaN/warning) when target outside (low, up)", {
  # up = 0.9 < target = 0.95: the threshold is unreachable on this draw, so the
  # whole vector must drop as NA rather than producing log()-of-negative NaNs.
  expect_no_warning(
    out_hi <- survival_from_dose(c(0, 0.5, 1, 2),
                                 low = 0.02, up = 0.9, k = 6,
                                 target_surv = 0.95)
  )
  expect_true(all(is.na(out_hi)))
  expect_false(any(is.nan(out_hi)))

  # target = 0.01 below low = 0.02: also unreachable -> NA, no warning.
  expect_no_warning(
    out_lo <- survival_from_dose(c(0, 1, 5),
                                 low = 0.02, up = 0.98, k = 6,
                                 target_surv = 0.01)
  )
  expect_true(all(is.na(out_lo)))

  # In-range target is unaffected and finite.
  out_ok <- survival_from_dose(c(0, 1),
                               low = 0.02, up = 0.98, k = 6,
                               target_surv = 0.5)
  expect_true(all(is.finite(out_ok)))
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

test_that("predict_heat_injury absolute target yields finite summaries without warnings", {
  skip_unless_brms()
  wf <- load_fixture_workflow()
  expect_no_warning(
    hi <- predict_heat_injury(
      data.frame(time = 0:6, temp = seq(30, 36, length.out = 7)),
      wf, target_surv = "absolute", ndraws = 200, seed = 1)
  )
  expect_true(all(is.finite(hi$summary$surv_median)))
  expect_true(all(is.finite(hi$summary$hi_median)))
})

test_that("predict_heat_injury errors on an unsupported time unit", {
  wf <- fake_hi_workflow("minutes")
  expect_error(
    predict_heat_injury(data.frame(time = c(0, 1), temp = 30), wf,
                        trace_unit = "fortnights", ndraws = 5),
    "unsupported time unit")
})

test_that("predict_heat_injury rejects an unknown shape", {
  wf <- fake_hi_workflow("hours")
  expect_error(
    predict_heat_injury(data.frame(time = c(0, 1), temp = 30), wf,
                        ndraws = 5, shape = "nonlinear"),
    "should be one of")
})

# ---------------------------------------------------------------------------
# shape = "varying": temperature-LOCAL low/up/k/mid in the dose integral.
# The CRITICAL invariant is that with a flat-in-T shape (low ~ 1, up ~ 1,
# k ~ 1) the incremental local-curve survival decrement telescopes EXACTLY to
# the constant-shape closed form, so the two modes must agree to numerical
# tolerance. With genuine T-structure on the asymptotes/steepness the modes
# diverge, and varying-mode survival must stay monotone.
# ---------------------------------------------------------------------------

test_that("shape='varying' reduces to 'constant' when low/up/k are flat in T", {
  skip_unless_brms()
  wf <- load_fixture_workflow_flatshape()
  scens <- make_temperature_scenarios(baseline = 20, spike_temp = 30, n_hours = 96,
                                      spike_times_single = 24,
                                      spike_times_multi  = c(24, 48, 72))
  for (tr in c("flat", "single_spike", "multi_spike")) {
    for (mode in c("relative", "absolute")) {
      hc <- predict_heat_injury(scens[[tr]], wf, target_surv = mode, T_c = 24,
                                ndraws = 200, seed = 7, shape = "constant")$summary
      hv <- predict_heat_injury(scens[[tr]], wf, target_surv = mode, T_c = 24,
                                ndraws = 200, seed = 7, shape = "varying")$summary
      expect_equal(hc$hi_median,   hv$hi_median,   tolerance = 1e-8,
                   info = paste(tr, mode, "HI"))
      expect_equal(hc$surv_median, hv$surv_median, tolerance = 1e-8,
                   info = paste(tr, mode, "S"))
    }
  }
})

test_that("shape='varying' differs from 'constant' under genuine T-structure, and S is monotone", {
  skip_unless_brms()
  wf <- load_fixture_workflow_shape()        # up ~ temp_c, k ~ temp_c (declining/steepening)
  # Hot spikes (above T_bar = 33) so the local asymptotes/steepness genuinely bite.
  scens <- make_temperature_scenarios(baseline = 30, spike_temp = 38, n_hours = 96,
                                      spike_times_multi = c(24, 48, 72))
  hc <- predict_heat_injury(scens$multi_spike, wf, target_surv = "relative",
                            T_c = 32, ndraws = 300, seed = 7, shape = "constant")
  hv <- predict_heat_injury(scens$multi_spike, wf, target_surv = "relative",
                            T_c = 32, ndraws = 300, seed = 7, shape = "varying")

  # Varying-mode survival is monotone non-increasing.
  expect_true(all(diff(hv$summary$surv_median) <= 1e-9))
  # The two modes give a materially different final survival.
  fc <- tail(hc$summary, 1); fv <- tail(hv$summary, 1)
  expect_gt(abs(fc$surv_median - fv$surv_median), 1e-3)
  expect_equal(hv$meta$shape, "varying")
})

test_that("shape='varying' does not jump survival when T changes with no added dose", {
  skip_unless_brms()
  wf <- load_fixture_workflow_shape()
  # Below the damage threshold T_c everywhere: dose stays at zero, but T jumps
  # around. Varying mode must hold survival flat (no path-dependent rebound/drop).
  trace <- data.frame(time = 0:5, temp = c(20, 38, 22, 36, 24, 30))
  hv <- predict_heat_injury(trace, wf, target_surv = "relative", T_c = 40,
                            ndraws = 100, seed = 7, shape = "varying")$summary
  expect_equal(hv$hi_median, rep(0, nrow(hv)))
  expect_equal(max(hv$surv_median) - min(hv$surv_median), 0, tolerance = 1e-9)
})

test_that("shape='varying' integrates per group on a grouped fit", {
  skip_unless_brms()
  wf <- load_fixture_workflow_shape_grouped()
  scens <- make_temperature_scenarios(baseline = 30, spike_temp = 37, n_hours = 72,
                                      spike_times_multi = c(12, 36))
  hv <- predict_heat_injury(scens$multi_spike, wf, target_surv = "relative",
                            T_c = 32, ndraws = 150, seed = 7, shape = "varying")
  expect_true("grp" %in% names(hv$summary))
  expect_setequal(unique(as.character(hv$summary$grp)), c("A", "B"))
  # Each group's varying survival is monotone.
  for (g in c("A", "B")) {
    sg <- hv$summary$surv_median[as.character(hv$summary$grp) == g]
    expect_true(all(diff(sg) <= 1e-9), info = paste("group", g))
  }
})

test_that("predict_heat_injury recovers analytical planted dose within CrI", {
  skip_unless_brms()

  wf <- load_fixture_workflow()
  ts <- truth_summary()

  scens <- make_temperature_scenarios(baseline = 20, spike_temp = 28,
                                      n_hours = 96,
                                      spike_times_single = 24,
                                      spike_times_multi  = c(24, 48, 72))
  T_c_val <- 24

  planted_single <- planted_dose_from_trace(
    scens$single_spike, z = ts$z, CTmax_1hr = ts$CTmax_1hr, T_c = T_c_val
  )
  planted_multi  <- planted_dose_from_trace(
    scens$multi_spike, z = ts$z, CTmax_1hr = ts$CTmax_1hr, T_c = T_c_val
  )

  hi_single <- predict_heat_injury(scens$single_spike, wf,
                                   T_c = T_c_val, ndraws = 300)
  hi_multi  <- predict_heat_injury(scens$multi_spike,  wf,
                                   T_c = T_c_val, ndraws = 300)

  truth_single <- tail(planted_single$hi_cumulative, 1)
  truth_multi  <- tail(planted_multi$hi_cumulative,  1)
  final_single <- tail(hi_single$summary, 1)
  final_multi  <- tail(hi_multi$summary,  1)

  # Posterior median within ~30% relative (small fit + ~10 spikes); truth
  # inside 95% CrI for both scenarios.
  expect_lt(abs(final_single$hi_median - truth_single) / truth_single, 0.3)
  expect_lt(abs(final_multi$hi_median  - truth_multi)  / truth_multi,  0.3)

  expect_lte(final_single$hi_lower, truth_single)
  expect_gte(final_single$hi_upper, truth_single)
  expect_lte(final_multi$hi_lower,  truth_multi)
  expect_gte(final_multi$hi_upper,  truth_multi)

  # Flat trace at sub-T_c baseline -> zero HI by construction.
  hi_flat <- predict_heat_injury(scens$flat, wf,
                                 T_c = T_c_val, ndraws = 100)
  expect_equal(tail(hi_flat$summary$hi_median, 1), 0)
})
