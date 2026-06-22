# Fast, brms-free tests for the z derivation. The z math now lives in the pure
# engine helper tls_local_z(), which consumes already-evaluated logLT matrices
# (posterior_linpred-sourced in production). We test it here on SYNTHETIC logLT
# matrices built from known coefficients, so the closed-form identities are
# checked deterministically without a Stan fit. End-to-end derive_z() on a real
# fit is covered gated in test-extract_tdt.R / test-direct-fixture.R.

# Closed-form log10 LT at temperature(s) for per-draw coefficient sets.
synth_logLT <- function(temps, df, Tbar, bnd, mode = "relative", p = 0.5) {
  mc  <- outer(rep(1, nrow(df)), temps - Tbar)               # [nd x nT]
  mid <- df$b_mid_Intercept + df$b_mid_temp_c * mc
  if (mode == "relative") return(mid)
  ell <- bnd$low_min + bnd$low_w * stats::plogis(df$b_lowraw_Intercept + df$b_lowraw_temp_c * mc)
  u   <- bnd$up_min  + bnd$up_w  * stats::plogis(df$b_upraw_Intercept  + df$b_upraw_temp_c  * mc)
  k   <- exp(df$b_logk_Intercept + df$b_logk_temp_c * mc)
  mid + log((u - p) / (p - ell)) / k
}
# z(T) via tls_local_z from logLT evaluated at temps -/+ h.
zfrom <- function(temps, df, Tbar, bnd, mode = "relative", p = 0.5, h = 1e-3) {
  tls_local_z(synth_logLT(temps + h, df, Tbar, bnd, mode, p),
              synth_logLT(temps - h, df, Tbar, bnd, mode, p), h, temps)
}

# Analytic local z(T) = -1 / (d/dT log10 LT_p) for a single coefficient set.
analytic_local_z <- function(T, co, Tbar, bnd, p = 0.5) {
  mc <- T - Tbar
  sl <- stats::plogis(co$a0 + co$a1 * mc); su <- stats::plogis(co$c0 + co$c1 * mc)
  ell <- bnd$low_min + bnd$low_w * sl; u <- bnd$up_min + bnd$up_w * su
  k   <- exp(co$g0 + co$g1 * mc)
  ellp <- bnd$low_w * sl * (1 - sl) * co$a1; up <- bnd$up_w * su * (1 - su) * co$c1
  g <- log((u - p) / (p - ell)); gp <- up / (u - p) + ellp / (p - ell)
  -1 / (co$b1 + (gp - g * co$g1) / k)
}

test_that("tls_local_z relative gives z = -1 / b_mid_temp_c exactly, per draw", {
  set.seed(1); nd <- 400
  b1  <- stats::rnorm(nd, -0.18, 0.02)
  bnd <- compute_4pl_bounds(0, 1)
  df  <- data.frame(
    b_mid_Intercept = stats::rnorm(nd, 1.7, 0.02), b_mid_temp_c = b1,
    b_lowraw_Intercept = stats::rnorm(nd, -3.0, 0.05), b_lowraw_temp_c = stats::rnorm(nd, 0.05, 0.01),
    b_upraw_Intercept = stats::rnorm(nd, 2.2, 0.05),  b_upraw_temp_c = stats::rnorm(nd, -0.25, 0.01),
    b_logk_Intercept = stats::rnorm(nd, 1.3, 0.05),   b_logk_temp_c = stats::rnorm(nd, 0.12, 0.01))
  temps <- c(28, 30, 32, 34, 36)
  zr <- zfrom(temps, df, Tbar = mean(temps), bnd, mode = "relative")
  # central difference of a LINEAR mid is exact -> z = -1/b1 to machine precision
  expect_equal(zr$draws$z, -1 / b1, tolerance = 1e-12)
  expect_equal(unique(round(zr$local_summary$z_median, 8)),
               round(stats::median(-1 / b1), 8))
})

test_that("tls_local_z absolute local z(T) matches the analytic derivative", {
  set.seed(2); nd <- 300
  bnd <- compute_4pl_bounds(0, 1)
  df  <- data.frame(
    b_mid_Intercept = stats::rnorm(nd, 1.7, 0.02), b_mid_temp_c = stats::rnorm(nd, -0.18, 0.02),
    b_lowraw_Intercept = stats::rnorm(nd, -3.0, 0.05), b_lowraw_temp_c = stats::rnorm(nd, 0.05, 0.01),
    b_upraw_Intercept = stats::rnorm(nd, 2.2, 0.05),  b_upraw_temp_c = stats::rnorm(nd, -0.25, 0.01),
    b_logk_Intercept = stats::rnorm(nd, 1.3, 0.05),   b_logk_temp_c = stats::rnorm(nd, 0.12, 0.01))
  temps <- c(28, 32, 36); Tbar <- mean(temps)
  za <- zfrom(temps, df, Tbar, bnd, mode = "absolute")
  co_of <- function(i) list(b1 = df$b_mid_temp_c[i],
    a0 = df$b_lowraw_Intercept[i], a1 = df$b_lowraw_temp_c[i],
    c0 = df$b_upraw_Intercept[i], c1 = df$b_upraw_temp_c[i],
    g0 = df$b_logk_Intercept[i], g1 = df$b_logk_temp_c[i])
  for (j in seq_along(temps)) {
    an <- vapply(seq_len(nd), function(i) analytic_local_z(temps[j], co_of(i), Tbar, bnd), numeric(1))
    fd <- za$local_draws$z[za$local_draws$temp == temps[j]]
    expect_equal(fd, an, tolerance = 1e-5)             # central diff matches to O(h^2)
  }
})

test_that("tls_local_z absolute reduces to -1/b_mid_temp_c when shape is constant in T", {
  set.seed(3); nd <- 200
  b1  <- stats::rnorm(nd, -0.18, 0.02)
  bnd <- compute_4pl_bounds(0, 1)
  df  <- data.frame(
    b_mid_Intercept = stats::rnorm(nd, 1.7, 0.02), b_mid_temp_c = b1,
    b_lowraw_Intercept = stats::rnorm(nd, -3.0, 0.05), b_lowraw_temp_c = 0,
    b_upraw_Intercept = stats::rnorm(nd, 2.2, 0.05),  b_upraw_temp_c = 0,
    b_logk_Intercept = stats::rnorm(nd, 1.3, 0.05),   b_logk_temp_c = 0)
  temps <- c(28, 30, 32, 34)
  za <- zfrom(temps, df, Tbar = mean(temps), bnd, mode = "absolute")
  expect_equal(za$draws$z, -1 / b1, tolerance = 1e-5)  # flat correction -> slope = b1
})

test_that("derive_z errors on an unfitted workflow", {
  wf <- structure(list(fit = NULL), class = "bayes_tls")
  expect_error(derive_z(wf), "Fit the model first")
})
