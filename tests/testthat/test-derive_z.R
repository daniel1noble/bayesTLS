# Fast, brms-free unit tests for derive_z(). derive_z() only needs
# posterior::as_draws_df(workflow$fit), so we hand it a fake `bayes_tls`
# workflow whose `$fit` is a draws_df of KNOWN coefficients. This lets us check
# the posterior z calculation against closed-form truth deterministically,
# without fitting a Stan model. The brms-gated recovery test lives in
# test-extract_tdt.R.

# Closed-form 4PL parametrisation used by the package (disjoint bounds):
#   ell(T) = low_min + low_w * plogis(lowraw(T))
#   u(T)   = up_min  + up_w  * plogis(upraw(T))
#   k(T)   = exp(logk(T)),   mid(T) = b0 + b1 (T - Tbar)
# Absolute-threshold p: log10 LT_p(T) = mid + (1/k) log((u - p)/(p - ell)).
# Analytic local z(T) = -1 / (d/dT log10 LT_p) for a single coefficient set.
analytic_local_z <- function(T, co, Tbar, bnd, p = 0.5) {
  mc <- T - Tbar
  sl <- stats::plogis(co$a0 + co$a1 * mc)
  su <- stats::plogis(co$c0 + co$c1 * mc)
  ell <- bnd$low_min + bnd$low_w * sl
  u   <- bnd$up_min  + bnd$up_w  * su
  k   <- exp(co$g0 + co$g1 * mc)
  ellp <- bnd$low_w * sl * (1 - sl) * co$a1
  up   <- bnd$up_w  * su * (1 - su) * co$c1
  g    <- log((u - p) / (p - ell))
  gp   <- up / (u - p) + ellp / (p - ell)
  slope <- co$b1 + (gp - g * co$g1) / k
  -1 / slope
}

# Build a fake fitted bayes_tls whose $fit is a draws_df of known coefficients.
fake_workflow <- function(draws_df, temps, Tbar, bnd) {
  wf <- list(
    fit  = posterior::as_draws_df(draws_df),
    meta = list(temp_mean = Tbar, bounds = bnd),
    data = data.frame(temp = temps)
  )
  class(wf) <- "bayes_tls"
  wf
}

test_that("derive_z relative gives z = -1 / b_mid_temp_c exactly, per draw", {
  set.seed(1)
  nd  <- 400
  b1  <- stats::rnorm(nd, -0.18, 0.02)
  bnd <- compute_4pl_bounds(0, 1)
  df  <- data.frame(
    b_mid_Intercept    = stats::rnorm(nd, 1.7, 0.02),
    b_mid_temp_c       = b1,
    b_lowraw_Intercept = stats::rnorm(nd, -3.0, 0.05),
    b_lowraw_temp_c    = stats::rnorm(nd, 0.05, 0.01),
    b_upraw_Intercept  = stats::rnorm(nd, 2.2, 0.05),
    b_upraw_temp_c     = stats::rnorm(nd, -0.25, 0.01),
    b_logk_Intercept   = stats::rnorm(nd, 1.3, 0.05),
    b_logk_temp_c      = stats::rnorm(nd, 0.12, 0.01)
  )
  temps <- c(28, 30, 32, 34, 36)
  wf <- fake_workflow(df, temps, Tbar = mean(temps), bnd = bnd)

  zr <- derive_z(wf, target_surv = "relative")
  # Per-draw identity (sort both: derive_z filters to finite, order preserved).
  expect_equal(zr$draws$z, -1 / b1, tolerance = 1e-12)
  # Relative z is constant in temperature: every local z equals -1/b1.
  expect_equal(unique(round(zr$local_summary$z_median, 8)),
               round(stats::median(-1 / b1), 8))
})

test_that("derive_z absolute local z(T) matches the analytic derivative", {
  set.seed(2)
  nd  <- 300
  bnd <- compute_4pl_bounds(0, 1)
  df  <- data.frame(
    b_mid_Intercept    = stats::rnorm(nd, 1.7, 0.02),
    b_mid_temp_c       = stats::rnorm(nd, -0.18, 0.02),
    b_lowraw_Intercept = stats::rnorm(nd, -3.0, 0.05),
    b_lowraw_temp_c    = stats::rnorm(nd, 0.05, 0.01),
    b_upraw_Intercept  = stats::rnorm(nd, 2.2, 0.05),
    b_upraw_temp_c     = stats::rnorm(nd, -0.25, 0.01),
    b_logk_Intercept   = stats::rnorm(nd, 1.3, 0.05),
    b_logk_temp_c      = stats::rnorm(nd, 0.12, 0.01)
  )
  temps <- c(28, 32, 36)
  Tbar  <- mean(temps)
  wf    <- fake_workflow(df, temps, Tbar = Tbar, bnd = bnd)

  za <- derive_z(wf, target_surv = "absolute", temp_grid = temps)

  # Analytic local z per (draw, temp), compared to derive_z's finite-difference.
  for (j in seq_along(temps)) {
    co_all <- lapply(seq_len(nd), function(i) list(
      b0 = df$b_mid_Intercept[i], b1 = df$b_mid_temp_c[i],
      a0 = df$b_lowraw_Intercept[i], a1 = df$b_lowraw_temp_c[i],
      c0 = df$b_upraw_Intercept[i], c1 = df$b_upraw_temp_c[i],
      g0 = df$b_logk_Intercept[i], g1 = df$b_logk_temp_c[i]
    ))
    an <- vapply(co_all, function(co)
      analytic_local_z(temps[j], co, Tbar, bnd), numeric(1))
    fd <- za$local_draws$z[za$local_draws$temp == temps[j]]
    # Central finite difference (h = 1e-3) matches analytic to O(h^2).
    expect_equal(fd, an, tolerance = 1e-5)
  }

  # Pooled z is the per-draw mean of the local z over temp_grid.
  pooled_manual <- rowMeans(sapply(temps, function(T)
    vapply(seq_len(nd), function(i) analytic_local_z(T, list(
      b0 = df$b_mid_Intercept[i], b1 = df$b_mid_temp_c[i],
      a0 = df$b_lowraw_Intercept[i], a1 = df$b_lowraw_temp_c[i],
      c0 = df$b_upraw_Intercept[i], c1 = df$b_upraw_temp_c[i],
      g0 = df$b_logk_Intercept[i], g1 = df$b_logk_temp_c[i]), Tbar, bnd),
      numeric(1))))
  expect_equal(stats::median(za$draws$z), stats::median(pooled_manual),
               tolerance = 1e-4)
})

test_that("derive_z absolute reduces to -1/b_mid_temp_c when shape is constant in T", {
  set.seed(3)
  nd  <- 200
  b1  <- stats::rnorm(nd, -0.18, 0.02)
  bnd <- compute_4pl_bounds(0, 1)
  # No temperature effect on lowraw / upraw / logk -> flat correction term.
  df  <- data.frame(
    b_mid_Intercept    = stats::rnorm(nd, 1.7, 0.02),
    b_mid_temp_c       = b1,
    b_lowraw_Intercept = stats::rnorm(nd, -3.0, 0.05),
    b_lowraw_temp_c    = 0,
    b_upraw_Intercept  = stats::rnorm(nd, 2.2, 0.05),
    b_upraw_temp_c     = 0,
    b_logk_Intercept   = stats::rnorm(nd, 1.3, 0.05),
    b_logk_temp_c      = 0
  )
  temps <- c(28, 30, 32, 34)
  wf <- fake_workflow(df, temps, Tbar = mean(temps), bnd = bnd)

  za <- derive_z(wf, target_surv = "absolute", temp_grid = temps)
  # With a flat correction the absolute-threshold slope is exactly b1.
  expect_equal(za$draws$z, -1 / b1, tolerance = 1e-5)
})

test_that("derive_z uses all draws by default and ndraws subsamples", {
  set.seed(4)
  nd  <- 100
  b1  <- stats::rnorm(nd, -0.2, 0.03)
  bnd <- compute_4pl_bounds(0, 1)
  df  <- data.frame(
    b_mid_Intercept = stats::rnorm(nd, 1.5, 0.02),
    b_mid_temp_c    = b1
  )
  temps <- c(30, 32, 34)
  wf <- fake_workflow(df, temps, Tbar = mean(temps), bnd = bnd)

  # Default: every draw is used (no subsampling), so z matches -1/b1 exactly.
  z_all <- derive_z(wf, target_surv = "relative")
  expect_equal(nrow(z_all$draws), nd)
  expect_equal(z_all$draws$z, -1 / b1, tolerance = 1e-12)

  # ndraws subsamples to the requested count.
  z_sub <- derive_z(wf, target_surv = "relative", ndraws = 40)
  expect_equal(nrow(z_sub$draws), 40L)
})

test_that("derive_z errors on an unfitted workflow", {
  wf <- structure(list(fit = NULL), class = "bayes_tls")
  expect_error(derive_z(wf), "Fit the model first")
})
