# Tier-1 fast, deterministic, brms-free tests for the direct-mode z/CTmax/C(T)
# derivation in extract_tdt(). These feed hand-built coefficient draws to the
# internal helpers (tdt_z_from_pars / tdt_ctmax_from_pars / tdt_loglt) with a
# `direct` spec, so they need no Stan fit and run on every CI pass. They guard
# the silent-NA regression (a direct fit has no b_mid_* -> the old code returned
# all NA) and the C(T) conversion arithmetic against an independent oracle.

make_direct_pars <- function(np = 300, tbar = 30, ctmaxdev = 2, seed = 1) {
  set.seed(seed)
  data.frame(
    b_lowraw_Intercept   = rnorm(np, -3.2,    0.05),
    "b_lowraw_temp_c"     = rnorm(np,  0.10,   0.02),
    b_upraw_Intercept    = rnorm(np,  3.2,    0.05),
    "b_upraw_temp_c"      = rnorm(np, -0.10,   0.02),
    b_logk_Intercept     = rnorm(np,  log(2), 0.05),
    "b_logk_temp_c"       = rnorm(np,  0.05,   0.02),
    b_CTmaxdev_Intercept = rnorm(np,  ctmaxdev, 0.05),
    b_logz_Intercept     = rnorm(np,  log(3), 0.03),
    check.names = FALSE
  )
}
DIRECT <- function(threshold = "relative")
  list(ctmaxdev = "b_CTmaxdev_Intercept", logz = "b_logz_Intercept",
       fit_threshold = threshold, log10_tref = 0)

test_that("direct relative z/CTmax are finite and recover exp(logz)/Tbar+CTmaxdev (silent-NA guard)", {
  tbar <- 30; pars <- make_direct_pars(tbar = tbar, ctmaxdev = 2)
  bnd  <- compute_4pl_bounds(0, 1)
  ts   <- resolve_target_surv("relative")
  tg   <- seq(tbar - 4, tbar + 8, length.out = 25)
  zo <- tdt_z_from_pars(pars, tbar, bnd, ts, tg, local = FALSE, direct = DIRECT())
  co <- tdt_ctmax_from_pars(pars, tbar, bnd, ts, exposure_model = 1, tg, direct = DIRECT())
  # the bug returned all-NA; assert BOTH are finite (CTmax was the under-tested one)
  expect_true(is.finite(zo$summary$z_median))
  expect_true(is.finite(co$summary$temp_median))
  # and they recover the coefficients: z = exp(logz), CTmax = Tbar + CTmaxdev
  expect_equal(zo$summary$z_median, median(exp(pars$b_logz_Intercept)), tolerance = 1e-3)
  expect_equal(co$summary$temp_median, tbar + median(pars$b_CTmaxdev_Intercept), tolerance = 1e-2)
})

test_that("absolute - relative logLT equals the analytic C(T) (independent oracle)", {
  tbar <- 30; pars <- make_direct_pars(tbar = tbar)
  bnd  <- compute_4pl_bounds(0, 1)
  T    <- c(28, 30, 32)
  mrel <- tdt_loglt(pars, tbar, bnd, resolve_target_surv("relative"), T, DIRECT())
  mabs <- tdt_loglt(pars, tbar, bnd, resolve_target_surv("absolute"), T, DIRECT())
  # C(T) computed by hand from the same draws, independent of tdt_loglt's internals
  lp <- function(par) outer(pars[[paste0("b_", par, "_temp_c")]], T - tbar) +
                      pars[[paste0("b_", par, "_Intercept")]]
  l <- bnd$low_min + bnd$low_w * plogis(lp("lowraw"))
  u <- bnd$up_min  + bnd$up_w  * plogis(lp("upraw"))
  k <- exp(lp("logk"))
  Chand <- log((u - 0.5) / (0.5 - l)) / k
  expect_equal(mabs - mrel, Chand, tolerance = 1e-9)
})

test_that("sublethal bounds (up < p) make the absolute correction NaN, not silently wrong", {
  # up bounded to [0.501, 0.999] by default; force the curve below p = 0.5 by
  # asking for an absolute p above the attainable upper asymptote.
  tbar <- 30; pars <- make_direct_pars(tbar = tbar)
  bnd  <- compute_4pl_bounds(0, 1)
  ts99 <- resolve_target_surv(0.9999)   # p above max attainable up -> (u - p) < 0
  out  <- suppressWarnings(tdt_loglt(pars, tbar, bnd, ts99, c(30, 31), DIRECT()))
  expect_true(all(is.na(out)))
})

test_that("threshold genuinely changes the quantity on an asymmetric curve", {
  tbar <- 30; pars <- make_direct_pars(tbar = tbar)   # up != 1 - low, shape varies with T
  bnd  <- compute_4pl_bounds(0, 1)
  tg   <- seq(tbar - 4, tbar + 8, length.out = 25)
  z_rel <- tdt_z_from_pars(pars, tbar, bnd, resolve_target_surv("relative"), tg,
                           local = FALSE, direct = DIRECT())$summary$z_median
  z_abs <- tdt_z_from_pars(pars, tbar, bnd, resolve_target_surv("absolute"), tg,
                           local = FALSE, direct = DIRECT())$summary$z_median
  expect_false(isTRUE(all.equal(z_rel, z_abs)))   # relative != absolute when shape varies
})

# --- helper-function coefficient reads on a (fake) direct fit -----------------
# A minimal bayes_tls whose $fit is a draws_df with CTmaxdev/logz coefficients
# (no b_mid_*), exercising the direct branch of the coefficient-read helpers
# without a Stan fit.
fake_direct_workflow <- function(np = 200, tbar = 30, threshold = "relative",
                                 seed = 2) {
  set.seed(seed)
  df <- data.frame(
    b_lowraw_Intercept   = rnorm(np, -3.2,    0.05),
    b_upraw_Intercept    = rnorm(np,  3.2,    0.05),
    b_logk_Intercept     = rnorm(np,  log(5), 0.05),
    b_CTmaxdev_Intercept = rnorm(np,  2,      0.05),
    b_logz_Intercept     = rnorm(np,  log(3), 0.03)
  )
  structure(list(fit  = posterior::as_draws_df(df),
                 meta = list(temp_mean = tbar, bounds = compute_4pl_bounds(0, 1),
                             parameterization = "direct", threshold = threshold,
                             log10_tref = 0, duration_unit = "hours")),
            class = "bayes_tls")
}

test_that("extract_4pl_pars reconstructs the linear midpoint on a direct fit (silent-empty guard)", {
  wf <- fake_direct_workflow(tbar = 30)
  p  <- extract_4pl_pars(wf)
  expect_gt(nrow(p), 0)   # direct fits returned 0 rows before the rewire
  d  <- as.data.frame(posterior::as_draws_df(wf$fit))
  expect_equal(median(p$mid_temp),
               median(-1 / exp(d$b_logz_Intercept)), tolerance = 1e-9)
  expect_equal(median(p$mid_int),
               median(d$b_CTmaxdev_Intercept / exp(d$b_logz_Intercept)),
               tolerance = 1e-9)
})

test_that("tdt_parameter_table reports CTmax/z directly on a direct fit", {
  wf  <- fake_direct_workflow(tbar = 30)
  tab <- tdt_parameter_table(wf)
  expect_true(all(is.finite(tab$median)))
  expect_true(any(grepl("CTmax", tab$parameter)))
  d  <- as.data.frame(posterior::as_draws_df(wf$fit))
  expect_equal(tab$median[grepl("CTmax", tab$parameter)],
               median(30 + d$b_CTmaxdev_Intercept), tolerance = 1e-6)
  expect_equal(tab$median[tab$parameter == "z (°C)"],
               median(exp(d$b_logz_Intercept)), tolerance = 1e-6)
})

test_that("every single-condition helper redirects a grouped direct fit to tls() (finding 5)", {
  # Treatment-coded grouped direct fit: b_CTmaxdev_Intercept (reference level)
  # PLUS a group offset. The original guard (does the Intercept exist?) passed
  # this and silently reported the reference level only; tdt_is_grouped() catches
  # it coding-independently (>1 CTmaxdev/logz fixed-effect column => grouped).
  set.seed(3); np <- 100
  df <- data.frame(
    b_lowraw_Intercept     = rnorm(np, -3.2,    0.05),
    b_upraw_Intercept      = rnorm(np,  3.2,    0.05),
    b_logk_Intercept       = rnorm(np,  log(5), 0.05),
    b_CTmaxdev_Intercept   = rnorm(np,  2,      0.05),
    b_CTmaxdev_life_stageB = rnorm(np,  5,      0.05),   # second level => grouped
    b_logz_Intercept       = rnorm(np,  log(3), 0.03)
  )
  wf <- structure(list(fit  = posterior::as_draws_df(df),
                       data = data.frame(temp = c(30, 33, 36)),
                       meta = list(temp_mean = 33, bounds = compute_4pl_bounds(0, 1),
                                   parameterization = "direct", threshold = "relative",
                                   log10_tref = 0, duration_unit = "hours",
                                   group_vars = "life_stage")),
                  class = "bayes_tls")
  expect_error(tdt_parameter_table(wf),                           "tls")
  expect_error(extract_4pl_pars(wf),                              "tls")
  expect_error(derive_z(wf, temp_grid = c(30, 36)),               "tls")
  expect_error(derive_temperature_for_duration(wf, 1, c(28, 40)), "tls")
  expect_error(extract_tdt(wf),                                   "tls")
})

test_that("absolute-FIT direct reconstruction is correct in every threshold cell (findings 1 & 4)", {
  tbar <- 30; pars <- make_direct_pars(tbar = tbar)
  bnd  <- compute_4pl_bounds(0, 1)
  T    <- c(28, 30, 32)
  dabs <- DIRECT("absolute")           # p_fit defaults to 0.5 inside tdt_loglt
  # independent oracle: backbone and C(T; p) from the same draws
  lp  <- function(par) outer(pars[[paste0("b_", par, "_temp_c")]], T - tbar) +
                       pars[[paste0("b_", par, "_Intercept")]]
  l <- bnd$low_min + bnd$low_w * plogis(lp("lowraw"))
  u <- bnd$up_min  + bnd$up_w  * plogis(lp("upraw"))
  k <- exp(lp("logk"))
  bb <- 0 - (outer(rep(1, nrow(pars)), T - tbar) - pars$b_CTmaxdev_Intercept) /
            exp(pars$b_logz_Intercept)
  Cof <- function(p) log((u - p) / (p - l)) / k
  # abs-fit + relative request: backbone - C0(0.5), and FINITE (was all-NA: finding 1)
  rel <- tdt_loglt(pars, tbar, bnd, resolve_target_surv("relative"), T, dabs)
  expect_true(all(is.finite(rel)))
  expect_equal(rel, bb - Cof(0.5), tolerance = 1e-9)
  # abs-fit + absolute@0.5: backbone (the case masked at p = 0.5)
  expect_equal(tdt_loglt(pars, tbar, bnd, resolve_target_surv("absolute"), T, dabs),
               bb, tolerance = 1e-9)
  # abs-fit + custom absolute p = 0.3: backbone - C0(0.5) + C(0.3) (was = bb: finding 4)
  expect_equal(tdt_loglt(pars, tbar, bnd, resolve_target_surv(0.3), T, dabs),
               bb - Cof(0.5) + Cof(0.3), tolerance = 1e-9)
})

test_that("fit_4pl maps non-canonical duration_unit aliases consistently (finding 6)", {
  base <- data.frame(logd = log10(rep(c(1, 2, 4, 8), 3)),
                     temp_c = rep(c(-2, 0, 2), each = 4), n_surv = 5L, n_total = 10L)
  l10 <- function(u) {
    d <- base; attr(d, "tdt_meta") <- list(temp_mean = 35, duration_unit = u,
                                           response_type = "count")
    fit_4pl(d, ctmax = ~ 1, fit = FALSE)$meta$log10_tref
  }
  expect_equal(l10("h"),     l10("hours"))     # alias == canonical (was tm=1 fallback)
  expect_equal(l10("hr"),    l10("hours"))
  expect_equal(l10("Hours"), l10("hours"))     # case-insensitive
  expect_equal(l10("min"),   l10("minutes"))
})

test_that("derive_z is wired for direct fits (finding 2)", {
  wf <- fake_direct_workflow(tbar = 30)
  d  <- as.data.frame(posterior::as_draws_df(wf$fit))
  tg <- seq(26, 38, by = 0.5)
  zr <- derive_z(wf, temp_grid = tg)                          # relative (default)
  za <- derive_z(wf, temp_grid = tg, target_surv = "absolute")
  expect_true(is.finite(zr$summary$z_median))                 # was all-NA
  expect_equal(zr$summary$z_median, median(exp(d$b_logz_Intercept)), tolerance = 1e-3)
  expect_true(is.finite(za$summary$z_median))
})

test_that("derive_temperature_for_duration is wired for direct fits (finding 3)", {
  wf <- fake_direct_workflow(tbar = 30)
  d  <- as.data.frame(posterior::as_draws_df(wf$fit))
  tg <- seq(26, 40, by = 0.05)
  # relative request crosses the linear backbone at target = log10(1) = 0 = log10_tref,
  # so CTmax = Tbar + CTmaxdev exactly; was Tbar + target/0 = Inf before the wiring.
  out <- derive_temperature_for_duration(wf, exposure_duration = 1, temp_grid = tg,
                                         target_surv = "relative")
  expect_true(is.finite(out$summary$temp_median))
  expect_equal(out$summary$temp_median,
               median(30 + d$b_CTmaxdev_Intercept), tolerance = 0.05)
})

test_that("predict_heat_injury on a direct fit matches the analytical dose oracle", {
  # Build a direct fit with LT50(temp_mean) = 1 h: CTmaxdev = 0 (CTmax = temp_mean)
  # and log10_tref = log10(unit/hour) so the reference dose is one model-unit hour.
  z <- 5; tbar <- 30; ndraw <- 30
  bnd <- list(low_min = 0, low_w = 0.1, up_min = 0.9, up_w = 0.1)  # low=0.02, up=0.98
  unit_per_hour <- 60   # minutes-fitted model
  fit <- posterior::as_draws_df(data.frame(
    b_lowraw_Intercept   = rep(stats::qlogis((0.02 - bnd$low_min) / bnd$low_w), ndraw),
    b_upraw_Intercept    = rep(stats::qlogis((0.98 - bnd$up_min)  / bnd$up_w),  ndraw),
    b_logk_Intercept     = rep(log(6), ndraw),
    b_CTmaxdev_Intercept = rep(0, ndraw),
    b_logz_Intercept     = rep(log(z), ndraw)
  ))
  wf <- structure(list(fit = fit,
                       meta = list(temp_mean = tbar, duration_unit = "minutes",
                                   bounds = bnd, parameterization = "direct",
                                   threshold = "relative",
                                   log10_tref = log10(unit_per_hour))),
                  class = "bayes_tls")
  trace <- data.frame(time = 0:11,
                      temp = c(20, 22, 26, 30, 33, 34, 33, 30, 26, 24, 22, 20))
  hi <- predict_heat_injury(trace, wf, trace_unit = "hours", T_c = 25,
                            ndraws = 20)$summary
  pl <- planted_dose_from_trace(trace, z = z, CTmax_1hr = tbar, T_c = 25)
  expect_equal(hi$hi_median, pl$hi_cumulative, tolerance = 1e-4)
})

# --- seed reproducibility (D-T2): same seed -> identical, different seed -> not -
test_that("seed makes derive_z reproducible across the draw subsample", {
  wf <- fake_direct_workflow(tbar = 30)          # 200 draws
  tg <- seq(26, 38, by = 0.5)
  a <- derive_z(wf, temp_grid = tg, ndraws = 40, seed = 42)$summary$z_median
  b <- derive_z(wf, temp_grid = tg, ndraws = 40, seed = 42)$summary$z_median
  cc <- derive_z(wf, temp_grid = tg, ndraws = 40, seed = 7)$summary$z_median
  expect_identical(a, b)                          # same seed -> identical subsample
  expect_false(isTRUE(all.equal(a, cc)))          # different seed -> different
})

test_that("seed makes predict_heat_injury reproducible across the draw subsample", {
  wf <- fake_direct_workflow(tbar = 30)          # 200 draws, relative threshold
  tr <- data.frame(time = 0:6, temp = c(26, 28, 30, 32, 30, 28, 26))
  s1 <- predict_heat_injury(tr, wf, ndraws = 30, seed = 1)$summary$surv_median
  s2 <- predict_heat_injury(tr, wf, ndraws = 30, seed = 1)$summary$surv_median
  s3 <- predict_heat_injury(tr, wf, ndraws = 30, seed = 2)$summary$surv_median
  expect_identical(s1, s2)
  expect_false(isTRUE(all.equal(s1, s3)))
})
