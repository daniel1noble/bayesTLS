# `local` gate in tls_local_z() (engine pure-math layer): skipping the
# per-temperature breakdown must not change the pooled z, and must only suppress
# local_draws/local_summary. We feed SYNTHETIC logLT matrices evaluated at
# grid -/+ h (no brms fit), which is exactly what derive_z()/extract_tdt() pass.

grid <- c(30, 33, 36, 39)
Tbar <- 34.5
h    <- 1e-3

# A bent logLT(T) per draw (quadratic term) so local z(T) genuinely varies in T.
make_logLT <- function(temps, np = 300, seed = 42) {
  set.seed(seed)                                  # same draws each call -> consistent
  a     <- stats::rnorm(np, 2.0, 0.10)
  slope <- stats::rnorm(np, 0.15, 0.01)           # z ~ 1/slope > 0
  curve <- stats::rnorm(np, 0.010, 0.005)         # bends the curve
  tc    <- temps - Tbar
  a - outer(slope, tc) + outer(curve, tc^2)       # [np x length(temps)]
}
lp <- make_logLT(grid + h)
lm <- make_logLT(grid - h)

test_that("pooled z is identical with local on vs off", {
  z_off <- tls_local_z(lp, lm, h, grid, local = FALSE)
  z_on  <- tls_local_z(lp, lm, h, grid, local = TRUE)
  expect_identical(z_off$draws,   z_on$draws)
  expect_identical(z_off$summary, z_on$summary)
})

test_that("local block is NULL only when local = FALSE", {
  z_off <- tls_local_z(lp, lm, h, grid, local = FALSE)
  z_on  <- tls_local_z(lp, lm, h, grid, local = TRUE)
  expect_null(z_off$local_draws)
  expect_null(z_off$local_summary)
  expect_s3_class(z_on$local_draws,   "tbl_df")
  expect_s3_class(z_on$local_summary, "tbl_df")
  expect_equal(nrow(z_on$local_summary), length(grid))   # one row per temperature
})

test_that("default is local = TRUE", {
  z_def <- tls_local_z(lp, lm, h, grid)            # no local arg
  expect_s3_class(z_def$local_summary, "tbl_df")
})

test_that("local z genuinely varies across temperature on a bent curve", {
  z_on <- tls_local_z(lp, lm, h, grid, local = TRUE)
  expect_gt(stats::sd(z_on$local_summary$z_median), 0)
})

test_that("local z of a LINEAR logLT is constant and equals -1/slope (relative-mode exactness)", {
  set.seed(7); np <- 200; b1 <- stats::rnorm(np, 0.18, 0.02)
  lin <- function(temps) 1.5 - outer(b1, temps - Tbar)   # linear -> central diff exact
  z <- tls_local_z(lin(grid + h), lin(grid - h), h, grid, local = TRUE)
  # central difference of a linear function is exact up to fp cancellation (~1e-11
  # per draw at h = 1e-3); the relative-z median parity vs the closed form is ~1e-13.
  expect_equal(z$draws$z, 1 / b1, tolerance = 1e-6)        # -1/(-b1) = 1/b1
  expect_lt(stats::sd(z$local_summary$z_median), 1e-6)     # constant in T
})
