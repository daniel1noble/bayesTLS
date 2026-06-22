# The vectorised CTmax crossing in tls_invert_logLT() (engine pure-math layer)
# must reproduce the original per-draw stats::approx() inverse interpolation
# exactly (machine precision), including non-monotone bent curves routed to the
# fallback. We feed SYNTHETIC logLT matrices directly (no brms fit, no
# coefficient parsing), which is precisely the matrix tls_invert_logLT consumes.

grid <- seq(28, 42, by = 0.05)

# Reference = the original per-row order + approx implementation.
ref_ctmax <- function(M, grid, target) {
  vapply(seq_len(nrow(M)), function(i) {
    y <- M[i, ]; ok <- is.finite(y)
    if (sum(ok) < 2L) return(NA_real_)
    o <- order(y[ok])
    suppressWarnings(stats::approx(y[ok][o], grid[ok][o], xout = target)$y)
  }, numeric(1))
}

expect_invert_equivalent <- function(M) {
  target <- stats::median(M, na.rm = TRUE)              # ensure curves cross it
  ref    <- ref_ctmax(M, grid, target)
  got    <- tls_invert_logLT(M, target, grid)           # per-draw Tc (NA where no crossing)
  fin    <- is.finite(ref)
  expect_equal(got[fin], ref[fin], tolerance = 1e-10)
  invisible(M)
}

test_that("monotone curves: vectorised tls_invert_logLT == approx-loop (machine precision)", {
  set.seed(1); np <- 500; gc <- grid - mean(grid)
  inter <- stats::rnorm(np, 1.5, 0.3)
  slope <- stats::rnorm(np, 0.18, 0.02)                 # > 0 -> logLT decreasing in T
  M <- inter - outer(slope, gc)                          # [np x nT], strictly monotone
  expect_invert_equivalent(M)
  dM <- M[, -1] - M[, -ncol(M)]
  expect_equal(sum(rowSums(dM >= 0) > 0), 0L)            # all monotone (vectorised path)
})

test_that("non-monotone bent curves route to the exact fallback and still match", {
  set.seed(2); np <- 400; gc <- grid - mean(grid)
  inter <- stats::rnorm(np, 1.5, 0.3)
  slope <- stats::rnorm(np, 0.05, 0.01)
  curve <- stats::rnorm(np, 0.03, 0.01)                  # quadratic term bends logLT
  M <- inter - outer(slope, gc) + outer(curve, gc^2)
  expect_invert_equivalent(M)
  dM <- M[, -1] - M[, -ncol(M)]
  expect_gt(sum(rowSums(dM >= 0) > 0), 0L)               # fallback actually exercised
})

test_that("rows with < 2 finite values yield NA (no crossing)", {
  M <- matrix(NA_real_, nrow = 3, ncol = length(grid))
  M[1, ] <- seq(2, -2, length.out = length(grid))        # crosses 0
  expect_true(is.finite(tls_invert_logLT(M, 0, grid)[1]))
  expect_true(all(is.na(tls_invert_logLT(M, 0, grid)[2:3])))
})

test_that("a flat (degenerate) draw yields NA without aborting the whole call", {
  # A flat draw has duplicate y; the old approx() fallback errored ("need at
  # least two non-NA values to interpolate") and took down the entire vectorised
  # inversion. It must now return NA for that row and finite for the others.
  M <- rbind(seq(2, -2, length.out = length(grid)),       # normal crossing
             rep(1.5, length(grid)))                       # flat -> no crossing
  got <- expect_no_error(tls_invert_logLT(M, target = 0, grid))
  expect_true(is.finite(got[1]))
  expect_true(is.na(got[2]))
})
