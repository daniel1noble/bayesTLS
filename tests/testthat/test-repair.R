# Tests for the Sharpe-Schoolfield repair kernel (R/repair.R). These are
# brms-free: repair_rate_schoolfield() is pure math, so every assertion checks
# an analytic property of the Sharpe-Schoolfield (1981) equation rather than a
# smoke "does it run" check.

# Parameters reused across tests: optimum near 17 degC.
rp <- list(TA = 14065, TAL = 50000, TAH = 120000,
           TL = 10.5 + 273.15, TH = 22.5 + 273.15, TREF = 17 + 273.15,
           r_ref = 0.01)

call_rate <- function(temp_celsius, pars = rp) {
  repair_rate_schoolfield(
    temp_celsius = temp_celsius,
    TA = pars$TA, TAL = pars$TAL, TAH = pars$TAH,
    TL = pars$TL, TH = pars$TH, TREF = pars$TREF, r_ref = pars$r_ref)
}

test_that("rate equals an independent evaluation of the Schoolfield equation", {
  # Independent reimplementation written from the published formula
  # (numerator = Arrhenius rate; denominator = 1 + low-arm + high-arm).
  ss <- function(Tc) {
    Tk  <- Tc + 273.15
    num <- exp(rp$TA * (1 / rp$TREF - 1 / Tk))
    den <- 1 + exp(rp$TAL * (1 / Tk - 1 / rp$TL)) +
               exp(rp$TAH * (1 / rp$TH - 1 / Tk))
    rp$r_ref * num / den
  }
  grid <- seq(2, 32, by = 0.5)
  expect_equal(call_rate(grid), vapply(grid, ss, numeric(1)),
               tolerance = 1e-12)
})

test_that("with both inactivation arms switched off the curve is pure Arrhenius", {
  # TAL = TAH = 0 -> denominator is exactly 3, numerator is Arrhenius, so the
  # rate is strictly increasing in temperature. This is an independent check
  # (different closed form) of the numerator.
  flat <- list(TA = 12000, TAL = 0, TAH = 0,
               TL = rp$TL, TH = rp$TH, TREF = rp$TREF, r_ref = 0.02)
  grid <- seq(0, 40, by = 1)
  got  <- call_rate(grid, flat)
  Tk   <- grid + 273.15
  expected <- flat$r_ref * exp(flat$TA * (1 / flat$TREF - 1 / Tk)) / 3
  expect_equal(got, expected, tolerance = 1e-12)
  expect_true(all(diff(got) > 0))                     # strictly increasing
})

test_that("the rate is a unimodal thermal-performance curve with peak between TL and TH", {
  grid  <- seq(0, 35, by = 0.1)
  rates <- call_rate(grid)
  peak_T <- grid[which.max(rates)]
  # Optimum must sit inside the low/high inactivation midpoints.
  expect_gt(peak_T, rp$TL - 273.15)
  expect_lt(peak_T, rp$TH - 273.15)
  # Strictly increasing up to the peak, strictly decreasing after it.
  below <- rates[grid <  peak_T]
  above <- rates[grid >  peak_T]
  expect_true(all(diff(below) > 0))
  expect_true(all(diff(above) < 0))
})

test_that("rate scales linearly with r_ref", {
  g  <- seq(5, 30, by = 2.5)
  r1 <- call_rate(g, modifyList(rp, list(r_ref = 0.01)))
  r2 <- call_rate(g, modifyList(rp, list(r_ref = 0.05)))
  expect_equal(r2, 5 * r1, tolerance = 1e-12)
})

test_that("the result is vectorised and matches element-wise scalar calls", {
  g   <- c(5, 12, 17, 22, 30)
  vec <- call_rate(g)
  expect_length(vec, length(g))
  scal <- vapply(g, function(t) call_rate(t), numeric(1))
  expect_equal(vec, scal, tolerance = 1e-12)
})

test_that("non-finite rates (overflow) and any negatives are coerced to zero", {
  # A huge TA at a high temperature overflows the Arrhenius numerator to Inf;
  # Inf / finite is non-finite and must be returned as 0, not NA/Inf.
  blow <- modifyList(rp, list(TA = 1e9))
  out  <- call_rate(c(17, 200), blow)
  expect_true(all(is.finite(out)))
  expect_equal(out[2], 0)
  # Across a wide grid the kernel is never negative.
  expect_true(all(call_rate(seq(-20, 60, by = 1)) >= 0))
})

test_that("rate at TREF is suppressed only by the (small) inactivation arms", {
  # At T = TREF the Arrhenius numerator is exactly 1, so rate = r_ref / denom
  # where denom = 1 + low-arm + high-arm. With TL/TH ~5 degC either side the
  # arms are small but non-zero, so the rate is just below r_ref.
  tref_c <- rp$TREF - 273.15
  rate_ref <- call_rate(tref_c)
  den <- 1 + exp(rp$TAL * (1 / rp$TREF - 1 / rp$TL)) +
             exp(rp$TAH * (1 / rp$TH - 1 / rp$TREF))
  expect_equal(rate_ref, rp$r_ref / den, tolerance = 1e-12)
  expect_lt(rate_ref, rp$r_ref)
  expect_gt(rate_ref, 0.9 * rp$r_ref)                 # arms are small here
})
