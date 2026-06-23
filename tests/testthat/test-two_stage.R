# Tests for the canonical classical two-stage TDT pipeline (R/two_stage.R).

# A clean synthetic assay where the truth is known: log10(LT50) declines
# linearly with temperature (slope -0.15 -> z = 1/0.15 ~ 6.67), asymptotes near
# 0/1 so the binomial Stage-1 midpoint coincides with log10(LT50).
make_two_stage_data <- function(seed = 99) {
  set.seed(seed)
  temps <- c(30, 32, 34, 36, 38)
  durs  <- c(1, 5, 15, 45, 135, 405)
  grid  <- expand.grid(temp = temps, dur = durs, rep = 1:5)
  mid   <- 1.5 - 0.15 * (grid$temp - 34)             # log10(LT50)(temp)
  p     <- 0.02 + 0.96 / (1 + exp(8 * (log10(grid$dur) - mid)))
  grid$tot  <- 20L
  grid$surv <- rbinom(nrow(grid), grid$tot, p)
  grid[, c("temp", "dur", "surv", "tot")]
}

test_that("ts_stage1 (binomial) recovers per-temperature LT50 and flags validity", {
  d  <- make_two_stage_data()
  s1 <- ts_stage1(d, "temp", "dur", "surv", "tot", family = "binomial")
  expect_setequal(s1$temp, c(30, 32, 34, 36, 38))
  expect_true(all(c("log10_lt50", "se_log10_lt50", "slope",
                    "finite_ok", "bracket_ok", "stage1_ok") %in% names(s1)))
  expect_true(all(s1$stage1_ok))                     # clean data: all valid
  expect_true(all(s1$slope < 0))                     # survival falls with time
  # midpoint at 34 C should sit near the simulated mid = 1.5
  expect_equal(s1$log10_lt50[s1$temp == 34], 1.5, tolerance = 0.15)
})

test_that("ts_stage2 recovers z, CTmax and T_crit", {
  d  <- make_two_stage_data()
  s2 <- ts_stage2(ts_stage1(d, "temp", "dur", "surv", "tot"), t_ref = 60)
  expect_equal(s2$summary$z, 1 / 0.15, tolerance = 0.6)        # ~6.67
  expect_true(is.finite(s2$summary$CTmax_1hr))
  expect_true(s2$summary$T_crit < s2$summary$CTmax_1hr)        # T_crit below CTmax
  expect_equal(s2$summary$n_stage1, 5L)
  expect_equal(s2$summary$n_excluded, 0L)
})

test_that("ts_ci delta gives wider t-quantile than Normal intervals", {
  d  <- make_two_stage_data()
  s2 <- ts_stage2(ts_stage1(d, "temp", "dur", "surv", "tot"))
  ci <- ts_ci(s2, method = "delta")
  expect_true(ci$z$lower < ci$z$point && ci$z$point < ci$z$upper)
  # t quantiles (few residual df) are wider than Normal
  expect_lt(ci$z$lower_t, ci$z$lower)
  expect_gt(ci$z$upper_t, ci$z$upper)
})

test_that("ts_ci mvn returns summary and curve bands", {
  d  <- make_two_stage_data()
  s2 <- ts_stage2(ts_stage1(d, "temp", "dur", "surv", "tot"))
  ci <- ts_ci(s2, method = "mvn", temp_grid = c(30, 34, 38), seed = 1)
  expect_true(all(c("z_lower", "z_upper", "CTmax_lower", "CTmax_upper",
                    "Tcrit_lower", "Tcrit_upper") %in% names(ci$summary_ci)))
  expect_equal(nrow(ci$curve_ci), 3L)
  expect_true(all(ci$curve_ci$duration_lower <= ci$curve_ci$duration_upper))
})

test_that("ts_curve returns a monotone-decreasing LT line", {
  d  <- make_two_stage_data()
  cu <- ts_curve(ts_stage2(ts_stage1(d, "temp", "dur", "surv", "tot")),
                 temp_grid = c(30, 34, 38))
  expect_equal(nrow(cu), 3L)
  expect_true(all(diff(cu$duration_median) < 0))     # shorter LT at hotter T
})

test_that("degenerate data yields no fit and NA summaries, not an error", {
  # Only two distinct durations per temperature -> Stage 1 cannot fit.
  d <- data.frame(temp = rep(c(30, 34, 38), each = 4),
                  dur  = rep(c(1, 400), times = 6),
                  surv = c(20, 0, 20, 1, 19, 0, 20, 0, 18, 2, 20, 1),
                  tot  = 20)
  s1 <- ts_stage1(d, "temp", "dur", "surv", "tot")
  expect_true(all(!s1$stage1_ok))
  s2 <- ts_stage2(s1)
  expect_null(s2$fit)
  expect_true(is.na(s2$summary$z))
  expect_silent(ts_ci(s2, method = "delta"))
  expect_silent(ts_ci(s2, method = "mvn", temp_grid = c(30, 34, 38)))
})

test_that("T_crit is invariant to t_ref (anchored at 1 h); CTmax tracks t_ref", {
  s1   <- ts_stage1(make_two_stage_data(), "temp", "dur", "surv", "tot")
  s60  <- ts_stage2(s1, t_ref = 60)$summary
  s120 <- ts_stage2(s1, t_ref = 120)$summary
  s240 <- ts_stage2(s1, t_ref = 240)$summary
  # rate-multiplier T_crit is anchored at the 1 h CTmax -> independent of the
  # reporting reference t_ref
  expect_equal(s60$T_crit, s120$T_crit, tolerance = 1e-8)
  expect_equal(s60$T_crit, s240$T_crit, tolerance = 1e-8)
  # the reported CTmax does shift with t_ref (a later reference -> lower temp)
  expect_lt(s120$CTmax_1hr, s60$CTmax_1hr)
  # ts_ci T_crit bounds are likewise t_ref-invariant
  ci60  <- ts_ci(ts_stage2(s1, t_ref = 60),  method = "mvn", seed = 1)$summary_ci
  ci120 <- ts_ci(ts_stage2(s1, t_ref = 120), method = "mvn", seed = 1)$summary_ci
  expect_equal(ci60$Tcrit_lower, ci120$Tcrit_lower, tolerance = 1e-6)
  expect_equal(ci60$Tcrit_upper, ci120$Tcrit_upper, tolerance = 1e-6)
})

test_that("ts_stage1 beta-binomial path works when glmmTMB is available", {
  skip_if_not_installed("glmmTMB")
  d  <- make_two_stage_data()
  s1 <- ts_stage1(d, "temp", "dur", "surv", "tot", family = "betabinomial")
  expect_equal(nrow(s1), 5L)
  expect_true("phi" %in% names(s1))                  # overdispersion reported
  expect_true(all(s1$stage1_ok))
})
