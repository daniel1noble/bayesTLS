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
# shape = "varying": temperature-LOCAL low/up/k/mid. Each individual of
# survival rank p accumulates injury against its own threshold time LT_p(T)
# and dies at 1; survival is the rank whose injury has reached 1. The CRITICAL
# invariant is that with a flat-in-T shape (low ~ 1, up ~ 1, k ~ 1) this is
# EXACTLY the constant-shape closed form, with or without repair, so the two
# modes must agree to numerical tolerance. With genuine T-structure on the
# asymptotes/steepness the modes diverge, and varying-mode survival must stay
# monotone.
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

test_that("shape='varying' with repair does not charge re-accrued injury twice", {
  skip_unless_brms()
  wf <- load_fixture_workflow_flatshape()
  # Each day 1.5 h at 30 C, then 22.5 h at 20 C (below T_c) where repair
  # removes most of the injury; the next day re-accrues it. The shape is flat
  # in temperature, so the two shapes must agree for every draw, with either
  # mortality setting and either kind of target.
  tr <- data.frame(time = seq(0, 96, by = 0.25))
  tr$temp <- ifelse(tr$time %% 24 < 1.5, 30, 20)
  rp <- list(TA = 14065, TAL = 50000, TAH = 100000, TL = 15 + 273.15,
             TH = 35 + 273.15, TREF = 25 + 273.15, r_ref = 0.1)
  for (meth in c("individual", "cumulative")) {
    for (irr in c(TRUE, FALSE)) {
      for (tg in list("relative", 0.5)) {
        args <- list(tr, wf, target_surv = tg, T_c = 24, repair = TRUE,
                     repair_pars = rp, irreversible_mortality = irr,
                     ndraws = 100, seed = 7, save_draws = TRUE)
        hc <- do.call(predict_heat_injury, c(args, shape = "constant"))$draws
        hv <- do.call(predict_heat_injury,
                      c(args, shape = "varying", varying_method = meth))$draws
        lab <- paste(meth, "| irreversible =", irr, "| target =", tg)
        expect_equal(hv$survival, hc$survival, tolerance = 1e-8, info = lab)
        expect_equal(hv$hi, hc$hi, tolerance = 1e-8, info = lab)
      }
    }
  }
})

test_that("shape='varying', varying_method='cumulative' follows its definition, draw by draw", {
  skip_unless_brms()
  wf <- load_fixture_workflow_shape()            # up and k change with temperature
  temps <- c(30, 31.5, 33, 34.5, 36)
  tr <- data.frame(time = seq(0, 72, by = 0.5))
  tr$temp <- temps[(seq_len(nrow(tr)) - 1) %/% 3 %% 5 + 1]
  rp <- list(TA = 14065, TAL = 50000, TAH = 100000, TL = 15 + 273.15,
             TH = 35 + 273.15, TREF = 25 + 273.15, r_ref = 0.02)
  rr <- repair_rate_schoolfield(tr$temp, TA = rp$TA, TAL = rp$TAL, TAH = rp$TAH,
                                TL = rp$TL, TH = rp$TH, TREF = rp$TREF,
                                r_ref = rp$r_ref)
  loc <- extract_4pl_pars(wf, temps = temps)
  for (rep_on in c(FALSE, TRUE)) {
    for (irr in c(TRUE, FALSE)) {
      pk <- predict_heat_injury(tr, wf, repair = rep_on, repair_pars = if (rep_on) rp,
                                irreversible_mortality = irr, shape = "varying",
                                varying_method = "cumulative", ndraws = 20, seed = 3,
                                save_draws = TRUE)$draws
      for (d in unique(pk$.draw)) {
        L <- loc[loc$.draw == d, ]
        P <- L[match(round(tr$temp, 6), round(L$temp, 6)), c("low", "up", "k", "mid")]
        want <- ref_cumulative(tr$time, tr$temp, P, rep_rate = if (rep_on) rr,
                               irreversible = irr)
        expect_equal(pk$survival[pk$.draw == d], want, tolerance = 1e-9,
                     info = paste("repair", rep_on, "irreversible", irr, "draw", d))
      }
    }
  }
})

test_that("shape='varying': the target sets heat injury only; a passed target is Inf, with a warning", {
  skip_unless_brms()
  wf <- load_fixture_workflow_shape()   # planted up(T) falls below 0.85 above ~35.5 C
  tm   <- seq(0, 48, by = 0.25)
  hot  <- data.frame(time = tm, temp = ifelse(tm %% 24 >= 12 & tm %% 24 < 15, 37, 33))
  base <- data.frame(time = tm, temp = 33)
  run  <- function(tr, tg) {
    predict_heat_injury(tr, wf, target_surv = tg, shape = "varying",
                        ndraws = 100, seed = 11, save_draws = TRUE)
  }
  rel <- run(hot, "relative")
  expect_warning(ab <- run(hot, 0.85), "already passed")
  # Survival does not depend on the target (no repair here).
  expect_equal(ab$draws$survival, rel$draws$survival, tolerance = 1e-12)
  # Heat injury is infinite on the draws whose 37 C curve starts below 0.85.
  up37 <- extract_4pl_pars(wf, temps = 37)
  end  <- ab$draws[ab$draws$time == max(tm), ]
  passed <- up37$up[match(end$.draw, up37$.draw)] <= 0.85
  expect_true(any(passed))
  expect_true(all(is.infinite(end$hi[passed])))
  expect_true(all(is.finite(end$hi[!passed])))
  # Hotter afternoons lower survival on every draw (before the fix they raised it).
  b_end <- run(base, 0.85)$draws
  b_end <- b_end[b_end$time == max(tm), ]
  expect_true(all(end$survival < b_end$survival[match(end$.draw, b_end$.draw)]))
  # A target passed at the very first step no longer drops the draw.
  first <- data.frame(time = tm, temp = ifelse(tm < 1, 37, 33))
  expect_warning(f <- run(first, 0.85), "already passed")
  expect_false(anyNA(f$draws$survival))
  expect_equal(f$meta$ndraws, 100)
})

test_that("shape='varying' keeps its floor, and ignores order and the start of the trace", {
  skip_unless_brms()
  wf <- load_fixture_workflow_shape()
  run_s <- function(tr) {
    d <- predict_heat_injury(tr, wf, shape = "varying", ndraws = 100, seed = 5,
                             save_draws = TRUE)$draws
    d <- d[d$time == max(tr$time), ]
    d[order(d$.draw), ]
  }
  # An hour at 36 C, then a long run at 31 C: survival stays above the floor.
  hs <- run_s(data.frame(time = 0:400, temp = c(36, rep(31, 400))))
  low31 <- extract_4pl_pars(wf, temps = 31)
  low31 <- low31$low[match(hs$.draw, low31$.draw)]
  expect_true(all(hs$survival >= low31 - 1e-12))
  expect_true(all(hs$survival > 0))
  # Without repair, the order of temperatures and where the trace starts in a
  # daily cycle do not matter (whole cycles, same draws).
  cyc <- function(sg) {
    tm <- seq(0, 96, by = 0.25)
    data.frame(time = tm, temp = round(33 + sg * 3 * cos(2 * pi * tm / 24), 4))
  }
  expect_equal(run_s(cyc(-1))$survival, run_s(cyc(1))$survival, tolerance = 1e-10)
  blk <- function(a, b) data.frame(time = (0:144) / 6, temp = c(rep(a, 72), rep(b, 73)))
  expect_equal(run_s(blk(36, 31))$survival, run_s(blk(31, 36))$survival, tolerance = 1e-10)
})

test_that("shape='varying' follows the per-individual rule, draw by draw", {
  skip_unless_brms()
  wf <- load_fixture_workflow_shape()
  temps <- c(30, 31.5, 33, 34.5, 36)
  tr <- data.frame(time = seq(0, 72, by = 0.5))
  tr$temp <- temps[(seq_len(nrow(tr)) - 1) %/% 3 %% 5 + 1]
  rp <- list(TA = 14065, TAL = 50000, TAH = 100000, TL = 15 + 273.15,
             TH = 35 + 273.15, TREF = 25 + 273.15, r_ref = 0.02)
  rr <- repair_rate_schoolfield(tr$temp, TA = rp$TA, TAL = rp$TAL, TAH = rp$TAH,
                                TL = rp$TL, TH = rp$TH, TREF = rp$TREF,
                                r_ref = rp$r_ref)
  loc <- extract_4pl_pars(wf, temps = temps)
  bar <- hi_pars(wf)
  # 0.85 is already passed at 36 C on some draws (heat injury Inf there)
  for (tg in list("relative", 0.7, 0.85)) {
    for (irr in c(TRUE, FALSE)) {
      for (scl in c(TRUE, FALSE)) {
        pk <- suppressWarnings(predict_heat_injury(
          tr, wf, target_surv = tg, repair = TRUE, repair_pars = rp,
          repair_scales_with_survival = scl, irreversible_mortality = irr,
          shape = "varying", ndraws = 20, seed = 3, save_draws = TRUE))$draws
        for (d in unique(pk$.draw)) {
          L <- loc[loc$.draw == d, ]
          P <- L[match(round(tr$temp, 6), round(L$temp, 6)), c("low", "up", "k", "mid")]
          p <- if (is.numeric(tg)) tg else NULL
          conv <- 1
          if (!is.null(p)) {
            B <- bar[bar$.draw == d, ]
            conv <- 10^(log((B$up - p) / (p - B$low)) / B$k)
          }
          want <- ref_individual(tr$time, tr$temp, P, p_target = p, conv = conv,
                                 rep_rate = rr, scale = scl, irreversible = irr)
          got  <- pk[pk$.draw == d, ]
          lab  <- paste("target", tg, "irreversible", irr, "scaled", scl, "draw", d)
          expect_equal(got$dose, want$dose, tolerance = 1e-10, info = lab)
          expect_equal(got$survival, want$survival, tolerance = 1e-10, info = lab)
        }
      }
    }
  }
  pk <- suppressWarnings(predict_heat_injury(tr, wf, target_surv = 0.85,
                                             shape = "varying", ndraws = 20, seed = 3,
                                             save_draws = TRUE))$draws
  expect_true(any(is.infinite(pk$hi)))
})

test_that("meta$ndraws counts the draws used, and a dropped draw warns", {
  skip_unless_brms()
  wf <- load_fixture_workflow_shape()
  # A target above up(T_bar) on about half the draws: those draws have no
  # constant-shape prediction.
  p <- round(stats::median(extract_4pl_pars(wf, temps = 33)$up), 3)
  expect_warning(
    o <- predict_heat_injury(data.frame(time = 0:24, temp = 34), wf,
                             target_surv = p, ndraws = 100, seed = 1,
                             save_draws = TRUE),
    "gave no prediction")
  used <- sum(tapply(is.finite(o$draws$survival), o$draws$.draw, all))
  expect_lt(used, 100)
  expect_equal(o$meta$ndraws, used)
  # a dropped draw is blanked in the heat-injury columns too, even when the
  # repair term does not involve survival
  expect_warning(
    o2 <- predict_heat_injury(data.frame(time = 0:24, temp = 34), wf,
                              target_surv = p, ndraws = 100, seed = 1,
                              repair_scales_with_survival = FALSE,
                              save_draws = TRUE),
    "gave no prediction")
  expect_true(all(is.na(o2$draws$hi[is.na(o2$draws$survival)])))
  expect_equal(o2$meta$ndraws, used)
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
