# Vectorised absolute-CTmax crossing in tdt_ctmax_from_pars() must reproduce the
# original per-draw stats::approx() inverse interpolation exactly (to machine
# precision), including for non-monotone bent curves routed to the fallback.
# Synthetic coefficient draws -> no brms fit needed (fast unit test).

bnd  <- list(low_min = 0.001, low_w = 0.498, up_min = 0.501, up_w = 0.498)
Tbar <- 34.5
grid <- seq(28, 42, by = 0.05)
ts   <- list(mode = "absolute", prob = 0.5, label = "p=0.500")

# Reference = the original implementation (per-row order + approx).
ref_ctmax <- function(M, grid, target) {
  vapply(seq_len(nrow(M)), function(i) {
    y <- M[i, ]; ok <- is.finite(y)
    if (sum(ok) < 2L) return(NA_real_)
    o <- order(y[ok])
    suppressWarnings(stats::approx(y[ok][o], grid[ok][o], xout = target)$y)
  }, numeric(1))
}

# Compare the full per-draw Tc vector reconstructed from the (NA-filtered) draws
# tibble returned by tdt_ctmax_from_pars() against the reference.
expect_ctmax_equivalent <- function(pars) {
  M       <- tdt_loglt(pars, Tbar, bnd, ts, grid)
  target  <- stats::median(M, na.rm = TRUE)        # ensure the curves cross it
  expo    <- 10^target                              # tdt_ctmax uses log10(expo) = target
  ref     <- ref_ctmax(M, grid, target)
  got     <- tdt_ctmax_from_pars(pars, Tbar, bnd, ts, expo, grid)$draws
  # Same draws are finite, and their values match the reference to machine eps.
  expect_identical(got$.draw, which(is.finite(ref)))
  expect_equal(got$temp, ref[got$.draw], tolerance = 1e-10)
  invisible(M)
}

test_that("monotone curves: vectorised CTmax == approx-loop reference (machine precision)", {
  set.seed(1)
  np <- 500
  pars <- data.frame(
    b_mid_Intercept    = stats::rnorm(np,  1.5,  0.05),
    b_mid_temp_c       = stats::rnorm(np, -0.15, 0.01),
    b_lowraw_Intercept = stats::rnorm(np, -2.0,  0.10),
    b_lowraw_temp_c    = stats::rnorm(np,  0.00, 0.02),
    b_upraw_Intercept  = stats::rnorm(np,  1.5,  0.10),
    b_upraw_temp_c     = stats::rnorm(np, -0.05, 0.01),
    b_logk_Intercept   = stats::rnorm(np, log(8), 0.10),
    b_logk_temp_c      = stats::rnorm(np,  0.02, 0.01)
  )
  M <- expect_ctmax_equivalent(pars)
  # Sanity: this set really is all monotone-decreasing (vectorised path only).
  dM <- M[, -1] - M[, -ncol(M)]
  expect_equal(sum(rowSums(dM >= 0) > 0), 0L)
})

test_that("non-monotone bent curves are routed to the exact fallback and still match", {
  # Steep positive temp_c slope on the upper asymptote makes the asymmetry
  # correction rise with T fast enough to bend log10 LT(T) non-monotone.
  set.seed(2)
  np <- 400
  pars <- data.frame(
    b_mid_Intercept    = stats::rnorm(np,  1.5,  0.05),
    b_mid_temp_c       = stats::rnorm(np, -0.02, 0.01),  # nearly flat midpoint
    b_lowraw_Intercept = stats::rnorm(np, -2.0,  0.10),
    b_lowraw_temp_c    = stats::rnorm(np,  0.00, 0.02),
    b_upraw_Intercept  = stats::rnorm(np,  0.0,  0.20),
    b_upraw_temp_c     = stats::rnorm(np,  1.0,  0.10),  # steep up rise -> bends LT
    b_logk_Intercept   = stats::rnorm(np, log(8), 0.10),
    b_logk_temp_c      = stats::rnorm(np,  0.05, 0.02)
  )
  M  <- expect_ctmax_equivalent(pars)
  # Confirm the fallback path was actually exercised (>=1 non-monotone row).
  dM <- M[, -1] - M[, -ncol(M)]
  expect_gt(sum(rowSums(dM >= 0) > 0), 0L)
})

test_that("relative CTmax (closed form) is unaffected by the absolute-path change", {
  set.seed(3); np <- 200
  pars <- data.frame(b_mid_Intercept = stats::rnorm(np, 1.5, 0.05),
                     b_mid_temp_c    = stats::rnorm(np, -0.15, 0.01))
  ts_rel <- list(mode = "relative", prob = NA_real_, label = "(low+up)/2")
  got <- tdt_ctmax_from_pars(pars, Tbar, bnd, ts_rel, 60, grid)$draws
  ref <- Tbar + (log10(60) - pars$b_mid_Intercept) / pars$b_mid_temp_c
  expect_equal(got$temp, ref[is.finite(ref)], tolerance = 1e-12)
})
