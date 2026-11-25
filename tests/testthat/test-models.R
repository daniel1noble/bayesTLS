# Alternative dose-response models (log-normal, Weibull type 1 / 2) through the
# model registry (R/models.R), the engine tdt_fit() and the exported wrappers.
# Fast tests need no Stan; the gated fits at the end simulate from each model
# and fit it with the matching function.

m_std <- function() {
  standardize_data(simulate_tdt(seed = 1), temp = "T", duration = "t_hours",
                   n_total = "n_trials", n_surv = "y_alive",
                   duration_unit = "hours")
}
m_std2 <- function() {
  standardize_data(simulate_tdt_2group(seed = 2), temp = "T",
                   duration = "t_hours", n_total = "n_trials",
                   n_surv = "y_alive", duration_unit = "hours")
}
# R-side evaluation of a Stan-syntax expression from the registry.
eval_stan <- function(expr, vars) {
  env <- list2env(c(vars, list(Phi = stats::pnorm, inv_Phi = stats::qnorm,
                               inv_logit = stats::plogis, log2 = log2)),
                  parent = baseenv())
  eval(parse(text = expr), env)
}
MODELS <- c("log-logistic", "log-normal", "weibull1", "weibull2")
grid <- expand.grid(low = c(0, 0.1, 0.3), up = c(1, 0.9, 0.6),
                    k = c(0.7, 2, 5), p = c(0.3, 0.5, 0.7))
grid <- grid[grid$low < grid$p & grid$p < grid$up, ]

# ---- registry mathematics ----------------------------------------------------

test_that("every model is at the midpoint of the range at u = 0 and decreasing", {
  for (m in MODELS) {
    G <- tdt_model_spec(m)$G
    expect_equal(G(0), 0.5, tolerance = 1e-12, info = m)
    u <- seq(-4, 4, by = 0.25)          # strictly decreasing where not saturated
    expect_true(all(diff(G(u)) < 0), info = m)
    expect_true(all(diff(G(seq(-8, 8, by = 0.25))) <= 0), info = m)
    expect_lt(G(8), 1e-3); expect_gt(G(-8), 1 - 1e-3)
  }
  w1 <- tdt_model_spec("weibull1")$G; w2 <- tdt_model_spec("weibull2")$G
  u  <- seq(-4, 4, by = 0.5)
  expect_equal(w2(u), 1 - w1(-u), tolerance = 1e-12)        # mirror image
})

test_that("the Stan formula strings are the registry curves", {
  for (m in MODELS) {
    spec <- tdt_model_spec(m)
    for (i in seq_len(nrow(grid))) {
      low <- grid$low[i]; up <- grid$up[i]; k <- grid$k[i]
      logd <- c(0.5, 1.5, 2.5); mid <- 1.7
      stan <- eval_stan(spec$rhs("low", "up"),
                        list(low = low, up = up, logk = log(k), logd = logd,
                             mid = mid))
      expect_equal(stan, low + (up - low) * spec$G(k * (logd - mid)),
                   tolerance = 1e-12, info = m)
    }
  }
})

test_that("absolute corrections solve p = target for every model, in R and in Stan syntax", {
  for (m in MODELS) {
    spec <- tdt_model_spec(m)
    for (i in seq_len(nrow(grid))) {
      low <- grid$low[i]; up <- grid$up[i]; k <- grid$k[i]; p <- grid$p[i]
      mid <- 1.7
      corr <- spec$abs_corr(low, up, k, p)
      curve <- function(x) low + (up - low) * spec$G(k * (x - mid))
      expect_equal(curve(mid + corr), p, tolerance = 1e-10, info = m)
      stan_corr <- eval_stan(spec$abs_corr_stan(p),
                             list(low = low, up = up, logk = log(k)))
      expect_equal(stan_corr, corr, tolerance = 1e-10, info = m)
    }
    # symmetric asymptotes about the target: no correction
    expect_equal(spec$abs_corr(0.1, 0.9, 2, 0.5), 0, tolerance = 1e-12, info = m)
    expect_equal(spec$abs_corr(0, 1, 3, 0.5), 0, tolerance = 1e-12, info = m)
  }
})

test_that("slope factors, inflection offsets and hazards match the curves", {
  h <- 1e-5
  for (m in MODELS) {
    spec <- tdt_model_spec(m)
    G <- spec$G
    expect_equal((G(h) - G(-h)) / (2 * h), -spec$slope_factor,
                 tolerance = 1e-6, info = m)
    # the standard location is the inflection point: G'' = 0 at u = theta_offset
    u0 <- spec$theta_offset
    expect_lt(abs((G(u0 + 1e-3) - 2 * G(u0) + G(u0 - 1e-3)) / 1e-6), 1e-3)
    # hazard for low = 0, up = 1 equals -d log p / dt numerically
    k <- 3; mid <- log10(60)
    p_fun <- function(t) G(k * (log10(t) - mid))
    for (t in c(10, 30, 60, 120, 400)) {
      num <- -(log(p_fun(t * (1 + h))) - log(p_fun(t * (1 - h)))) / (2 * t * h)
      expect_equal(spec$hazard(t, p_fun(t), k, mid), num,
                   tolerance = 1e-6, info = m)
    }
  }
})

# ---- names, priors, builders ---------------------------------------------------

test_that("model names and aliases resolve; unknown names error", {
  expect_equal(tdt_match_model(NULL), "log-logistic")
  expect_equal(tdt_match_model("Log-Logistic"), "log-logistic")
  expect_equal(tdt_match_model("4pl"), "log-logistic")
  expect_equal(tdt_match_model("probit"), "log-normal")
  expect_equal(tdt_match_model("lognormal"), "log-normal")
  expect_equal(tdt_match_model("gompertz"), "weibull1")
  expect_equal(tdt_match_model("weibull"), "weibull1")
  expect_equal(tdt_match_model("Weibull-2"), "weibull2")
  expect_equal(tdt_match_model("frechet"), "weibull2")
  expect_error(tdt_match_model("cauchy"), "Unknown model")
  expect_error(tdt_match_model(c("a", "b")), "`model` must be one of")
})

test_that("logk prior means equalise the midpoint-slope prior across models", {
  expect_identical(tdt_logk_prior_mean("log-logistic"), log(2))
  for (m in MODELS) {
    f <- tdt_model_spec(m)$slope_factor
    expect_equal(exp(tdt_logk_prior_mean(m)) * f, 2 * 0.25, tolerance = 1e-12,
                 info = m)
  }
})

test_that("make_4pl_* are the log-logistic case of make_tdt_*", {
  std <- m_std()
  expect_identical(utils::capture.output(print(make_4pl_formula())),
                   utils::capture.output(print(make_tdt_formula(model = "log-logistic"))))
  expect_identical(as.data.frame(make_4pl_priors(std)),
                   as.data.frame(make_tdt_priors(model = "log-logistic", data = std)))
  expect_identical(utils::capture.output(print(make_4pl_formula(ctmax = ~ 1, z = ~ 1))),
                   utils::capture.output(print(make_tdt_formula(model = "log-logistic",
                                                                ctmax = ~ 1, z = ~ 1))))
})

test_that("other models change only the logk Intercept prior", {
  std <- m_std()
  base <- as.data.frame(make_4pl_priors(std))
  drop_logk_int <- function(d) {
    d <- d[!(d$nlpar == "logk" & d$coef == "Intercept"), ]; rownames(d) <- NULL; d
  }
  for (m in setdiff(MODELS, "log-logistic")) {
    pr <- as.data.frame(make_tdt_priors(model = m, data = std))
    expect_identical(drop_logk_int(pr), drop_logk_int(base), info = m)
    lk <- pr[pr$nlpar == "logk" & pr$coef == "Intercept", "prior"]
    expect_equal(lk, sprintf("normal(%.6f, 1)", tdt_logk_prior_mean(m)), info = m)
    # direct parameterisation too
    prd <- as.data.frame(make_tdt_priors(model = m, data = std, ctmax = ~ 1, z = ~ 1))
    expect_equal(prd[prd$nlpar == "logk" & prd$coef == "Intercept", "prior"],
                 sprintf("normal(%.6f, 1)", tdt_logk_prior_mean(m)), info = m)
  }
})

# ---- wrappers and metadata -------------------------------------------------------

test_that("fit_4pl / fit_lognormal / fit_weibull are the engine with the model set", {
  std <- m_std()
  same <- function(a, b) {
    expect_identical(utils::capture.output(print(a$formula)),
                     utils::capture.output(print(b$formula)))
    expect_identical(as.data.frame(a$prior), as.data.frame(b$prior))
    expect_identical(a$meta, b$meta)
  }
  same(fit_4pl(std, fit = FALSE), tdt_fit(model = "log-logistic", data = std, fit = FALSE))
  same(fit_lognormal(std, fit = FALSE), tdt_fit(model = "log-normal", data = std, fit = FALSE))
  same(fit_weibull(std, fit = FALSE), tdt_fit(model = "weibull1", data = std, fit = FALSE))
  same(fit_weibull(std, type = 2, fit = FALSE), tdt_fit(model = "weibull2", data = std, fit = FALSE))
  expect_equal(fit_4pl(std, fit = FALSE)$meta$model, "log-logistic")
  expect_equal(fit_lognormal(std, fit = FALSE)$meta$model, "log-normal")
  expect_equal(fit_weibull(std, fit = FALSE)$meta$model, "weibull1")
  expect_equal(fit_weibull(std, type = 2, fit = FALSE)$meta$model, "weibull2")
  expect_error(fit_weibull(std, type = 3, fit = FALSE), "`type` must be 1")
  # arguments forward: a direct, fixed-asymptote Weibull spec
  wf <- fit_weibull(m_std2(), ctmax = ~ 0 + grp, z = ~ 0 + grp,
                    fix_asymptotes = c(up = 1), fit = FALSE)
  expect_equal(wf$meta$parameterization, "direct")
  expect_equal(names(wf$meta$fix_asymptotes), "up")
  expect_output(print(wf), "Model:.*Weibull type 1")
})

test_that("make_tdt_priors() takes the data first, like make_4pl_priors()", {
  std <- m_std()
  expect_identical(as.data.frame(make_tdt_priors(std)),
                   as.data.frame(make_4pl_priors(std)))
  expect_identical(as.data.frame(make_tdt_priors(std, "weibull1")),
                   as.data.frame(make_tdt_priors(std, model = "weibull1")))
})

test_that("each model is recognised from its own formula and no other", {
  std <- m_std(); std2 <- m_std2()
  for (m in MODELS) {
    specs <- list(
      tdt_fit(model = m, data = std, fit = FALSE),
      tdt_fit(model = m, data = std, ctmax = ~ 1, z = ~ 1, fit = FALSE),
      tdt_fit(model = m, data = std, ctmax = ~ 1, z = ~ 1,
              threshold = "absolute", fit = FALSE),
      tdt_fit(model = m, data = std2, by = "grp", fit = FALSE),
      tdt_fit(model = m, data = std, fix_asymptotes = TRUE, fit = FALSE))
    for (wf in specs) {
      expect_equal(tdt_model_from_formula(wf$formula), m, info = m)
      # a bare brmsfit carries only its formula: recognised the same way
      fake <- structure(list(formula = wf$formula), class = "brmsfit")
      expect_equal(tdt_model(fake), m, info = m)
      # a workflow saved before models existed (no meta$model)
      old <- wf; old$meta$model <- NULL
      expect_equal(tdt_model(old), m, info = m)
    }
  }
  expect_true(is.na(tdt_model_from_formula(y ~ a + b * x)))
  expect_true(is.na(tdt_model_from_formula("not a formula")))
  expect_equal(tdt_model(structure(list(formula = y ~ a), class = "brmsfit")),
               "log-logistic")                        # hand-written: assumed
})

test_that("tdt_model() canonicalises meta$model; print() survives a bad label", {
  wf <- fit_4pl(m_std(), fit = FALSE)
  wf$meta$model <- "4pl"
  expect_equal(tdt_model(wf), "log-logistic")
  wf$meta$model <- "gompertz"
  expect_equal(tdt_model(wf), "weibull1")
  wf$meta$model <- "bogus"
  expect_error(tdt_model(wf), "Unknown model")
  expect_output(print(wf), "Model:.*bogus \\(unrecognised\\)")
})

test_that("tdt_model() defaults pre-model objects to the log-logistic", {
  wf <- fit_4pl(m_std(), fit = FALSE)
  wf$meta$model <- NULL
  expect_equal(tdt_model(wf), "log-logistic")
  expect_equal(tdt_model(list()), "log-logistic")
})

test_that("Stan code generates for every model and configuration", {
  skip_if_not_installed("brms")
  std <- m_std(); std2 <- m_std2()
  for (m in MODELS) {
    for (wf in list(tdt_fit(model = m, data = std, fit = FALSE),
                    tdt_fit(model = m, data = std, ctmax = ~ 1, z = ~ 1,
                            threshold = "absolute", fit = FALSE),
                    tdt_fit(model = m, data = std2, by = "grp", fit = FALSE))) {
      code <- suppressWarnings(brms::make_stancode(wf$formula, data = wf$data,
                                                   prior = wf$prior))
      expect_type(code, "character")
      expect_true(nchar(code) > 1000, info = m)
    }
  }
})

# ---- absolute thresholds and heat injury for every model (Phase 2) -----------

test_that("curve() is low + (up - low) * G(u) over the whole range", {
  u <- seq(-12, 12, by = 0.05)
  for (m in MODELS) {
    spec <- tdt_model_spec(m)
    expect_equal(spec$curve(0.1, 0.9, u), 0.1 + 0.8 * spec$G(u),
                 tolerance = 1e-12, info = m)
  }
})

test_that("tdt_abs_offset() is the model's correction, NA where unreachable", {
  low <- matrix(c(0, 0.1, 0.6, 0.2), 2)     # cells [1,2] and [2,2] cannot
  up  <- matrix(c(1, 0.9, 0.95, 0.4), 2)    # reach p = 0.5 (p < low, p > up)
  k   <- matrix(c(2, 3, 4, 5), 2)
  reach <- low < 0.5 & 0.5 < up
  for (m in MODELS) {
    off <- expect_silent(tdt_abs_offset(m, low, up, k, 0.5))
    expect_equal(dim(off), dim(low))
    expect_true(all(is.na(off[!reach])))
    expect_equal(off[reach],
                 tdt_model_spec(m)$abs_corr(low[reach], up[reach], k[reach], 0.5))
  }
  # log-logistic: exactly the arithmetic of the pre-models engine
  ratio <- (up - 0.5) / (0.5 - low)
  ratio[!is.finite(ratio) | ratio <= 0] <- NA_real_
  expect_identical(tdt_abs_offset("log-logistic", low, up, k, 0.5), log(ratio) / k)
})

test_that("each model's curve passes through the target at the threshold time", {
  for (m in MODELS) {
    spec <- tdt_model_spec(m)
    for (i in seq_len(nrow(grid))) {
      g <- grid[i, ]
      tau <- time_to_surv_threshold_4pl(30, g$p, g$low, g$up, g$k,
                                        mid_int = 1.7, mid_temp = -0.2,
                                        temp_mean = 30, model = m)
      expect_equal(spec$curve(g$low, g$up, g$k * (log10(tau) - 1.7)), g$p,
                   tolerance = 1e-10, info = m)
    }
    expect_equal(time_to_surv_threshold_4pl(30, "relative", 0.1, 0.9, 3, 1.7,
                                            -0.2, 30, model = m), 10^1.7)
    expect_true(is.na(time_to_surv_threshold_4pl(30, 0.95, 0.1, 0.9, 3, 1.7,
                                                 -0.2, 30, model = m)))
  }
  # log-logistic: exactly the arithmetic of the pre-models code
  expect_identical(
    time_to_surv_threshold_4pl(c(28, 31), 0.4, 0.1, 0.9, 3, 1.7, -0.2, 30),
    10^((1.7 - 0.2 * (c(28, 31) - 30)) + log((0.9 - 0.4) / (0.4 - 0.1)) / 3))
})

test_that("one dose reaches the threshold on every model's survival curve", {
  doses <- c(1e-3, 0.1, 0.5, 1, 2, 10)
  for (m in MODELS) {
    for (i in seq_len(nrow(grid))) {
      g <- grid[i, ]
      expect_equal(survival_from_dose(1, g$low, g$up, g$k, target_surv = g$p,
                                      model = m), g$p, tolerance = 1e-10,
                   info = m)
      expect_equal(survival_from_dose(1, g$low, g$up, g$k, model = m),
                   (g$low + g$up) / 2, tolerance = 1e-12, info = m)
      # bounded and non-increasing everywhere (far tails reach an asymptote in
      # floating point, where low + (up - low) * 1 can round one unit past up);
      # strictly decreasing in the middle
      s <- survival_from_dose(doses, g$low, g$up, g$k, model = m)
      expect_true(all(diff(s) <= 1e-12) &&
                    all(s >= g$low - 1e-12 & s <= g$up + 1e-12), info = m)
      sc <- survival_from_dose(c(0.5, 1, 2), g$low, g$up, g$k, model = m)
      expect_true(all(diff(sc) < 0), info = m)
    }
    expect_true(all(is.na(survival_from_dose(doses, 0.1, 0.9, 3,
                                             target_surv = 0.95, model = m))))
  }
  # log-logistic: exactly the arithmetic of the pre-models code
  old_ll <- function(dose, low, up, k, p) {
    c_t <- log((up - p) / (p - low)) / k
    low + (up - low) / (1 + exp(k * (log10(pmax(dose, 1e-12)) + c_t)))
  }
  expect_identical(survival_from_dose(doses, 0.1, 0.9, 3, 0.4),
                   old_ll(doses, 0.1, 0.9, 3, 0.4))
})

# ---- wiring on a real (cached) fit, no Stan ------------------------------------
# The git-tracked Beta fixture has the nlpars every model uses (lowraw, upraw,
# logk, mid); swapping its formula for another model's makes brms evaluate that
# model on the same draws. fixtures/ is .Rbuildignore'd, so R CMD check skips.

beta_fixture_wf <- function(model = "log-logistic", ...) {
  path <- testthat::test_path("fixtures", "sim_fit_beta_small.rds")
  skip_if_not(file.exists(path), "Beta fixture not available")
  std <- standardize_data(simulate_tdt_beta(seed = 20260521), temp = "T",
                          duration = "t_hours", proportion = "y_prop",
                          duration_unit = "hours")
  fit <- readRDS(path)
  wf  <- tdt_fit(model = model, data = std, fit = FALSE, ...)
  fit$formula <- wf$formula          # same nlpars; brms now evaluates `model`
  wf$fit <- fit
  wf
}

test_that("tls() reads a workflow's own asymptote bounds", {
  skip_if_not_installed("brms")
  wf <- beta_fixture_wf()
  z  <- function(w, ...) tls(w, params = "z", target_surv = 0.3, ...)$summary$median
  # default bounds: reading them from the workflow changes nothing
  expect_identical(tls(wf, params = "z", target_surv = 0.3)$summary,
                   tls(wf, params = "z", target_surv = 0.3, lower = 0,
                       upper = 1)$summary)
  # a workflow fitted with other bounds is read with those bounds
  wf2 <- wf
  wf2$meta$bounds <- compute_4pl_bounds(0.2, 1)
  wf2$meta$lower  <- 0.2
  expect_equal(z(wf2), z(wf, lower = 0.2))
  expect_false(isTRUE(all.equal(z(wf2), z(wf))))
  # one explicit bound: the other comes from the workflow
  expect_equal(z(wf2, upper = 1), z(wf2))
})

test_that("absolute thresholds and heat injury are wired to each model (no Stan)", {
  skip_if_not_installed("brms")
  ids <- 1:60
  for (m in c("log-normal", "weibull1", "weibull2")) {
    wf   <- beta_fixture_wf(m)
    fit  <- get_brmsfit(wf)
    tbar <- wf$meta$temp_mean
    expect_identical(tdt_model(fit), m)
    temps <- c(30, 33, 36)
    nd <- tls_build_grid(fit$data, temp_grid = temps - tbar)
    sp <- tls_eval_subpars(fit, nd, wf$meta$bounds, draw_ids = ids,
                           mode = "relative")
    # the swap took effect: brms evaluates model m's curve
    pp1 <- brms::posterior_linpred(fit, newdata = new_tdt_grid(wf, temps, 1),
                                   transform = TRUE, re_formula = NA,
                                   draw_ids = ids)
    expect_equal(pp1, tdt_model_spec(m)$curve(sp$low, sp$up, -sp$k * sp$mid),
                 tolerance = 1e-12, ignore_attr = TRUE, info = m)
    # engine closed form vs numerical inversion of brms's predicted curve
    dgrid <- 10^seq(-3, 3, length.out = 1200)
    g  <- new_tdt_grid(wf, temps = temps, durations = dgrid)
    pp <- brms::posterior_linpred(fit, newdata = g, transform = TRUE,
                                  re_formula = NA, draw_ids = ids)
    for (p in c(0.3, 0.7)) {
      eng <- suppressWarnings(tls_eval_subpars(fit, nd, wf$meta$bounds,
                                               draw_ids = ids, mode = "absolute",
                                               p = p))$logLT
      for (j in seq_along(temps)) {
        num <- threshold_x_by_draw(pp[, g$temp == temps[j], drop = FALSE],
                                   dgrid, p)
        ok  <- is.finite(num) & is.finite(eng[, j])
        expect_gt(mean(ok), 0.9)
        expect_lt(max(abs(log10(num[ok]) - eng[ok, j])), 1e-3)
      }
    }
    # heat injury: absolute dose = sum(dt / LT_p) with LT_p from the engine,
    # for both integrators; survival on a constant T_bar trace is brms's curve
    trace <- data.frame(time = seq(0, 1, by = 0.25),
                        temp = c(31, 31, 35, 35, 35))
    tr_nd <- tls_build_grid(fit$data, temp_grid = trace$temp - tbar)
    flat  <- data.frame(time = seq(0, 2, by = 0.25), temp = tbar)
    cpp   <- brms::posterior_linpred(fit, newdata = new_tdt_grid(wf, tbar,
                                                                 flat$time[-1]),
                                     transform = TRUE, re_formula = NA)
    for (shp in c("constant", "varying")) {
      hi <- predict_heat_injury(trace, wf, target_surv = 0.3, ndraws = 20,
                                seed = 1, shape = shp, save_draws = TRUE)$draws
      d   <- sort(unique(hi$.draw))
      spa <- tls_eval_subpars(fit, tr_nd, wf$meta$bounds, draw_ids = d,
                              mode = "absolute", p = 0.3)
      logtau <- if (shp == "varying") spa$logLT else {
        sp0 <- tls_eval_subpars(fit, tls_build_grid(fit$data, temp_grid = 0),
                                wf$meta$bounds, draw_ids = d,
                                mode = "absolute", p = 0.3)
        spa$mid + as.vector(sp0$logLT - sp0$mid)
      }
      want <- t(apply(10^(-logtau[, -nrow(trace), drop = FALSE]), 1,
                      function(r) cumsum(c(0, r * diff(trace$time)))))
      got  <- t(sapply(d, function(i) hi$dose[hi$.draw == i]))
      expect_equal(got, want, tolerance = 1e-10, ignore_attr = TRUE,
                   info = paste(m, shp))
      s <- predict_heat_injury(flat, wf, target_surv = 0.3, ndraws = 20,
                               seed = 1, shape = shp, save_draws = TRUE)$draws
      for (i in unique(s$.draw))
        expect_equal(s$survival[s$.draw == i][-1], unname(cpp[i, ]),
                     tolerance = 1e-8, info = paste(m, shp))
    }
  }
})

test_that("tls() reads CTmax and z off each model's own threshold line (no Stan)", {
  skip_if_not_installed("brms")
  # References, per draw, that never use the models' corrections or tls()'s
  # handling of t_ref. Relative target: the midpoint line brms evaluates
  # (derive_tdt_curve()), solved in temperature. Absolute target: brms's
  # predicted survival solved for the threshold time on a dense duration grid,
  # then in temperature; its z is the local slope of that line averaged over
  # tls()'s default grid (11 temperatures across the data), as tls() averages
  # it. The absolute CTmax is checked at 60 and 30 min, so t_ref is not only
  # ever the fit's own. The fixture is in hours. These tests need the
  # git-tracked Beta fixture, which the package build leaves out
  # (tests/testthat/fixtures is in .Rbuildignore), so they run in
  # devtools::test(), not in R CMD check.
  temps <- seq(30, 36, by = 0.25)
  ids   <- 1:60
  zg    <- seq(30, 36, length.out = 11)
  d     <- 0.1
  t_ct  <- seq(30, 34.5, by = 0.25)
  t_abs <- sort(unique(c(t_ct, zg - d, zg + d)))
  dgrid <- 10^seq(log10(0.02), log10(20), length.out = 400)
  pick  <- function(dd, q, draws) {
    x <- dd[dd$quantity == q, ]
    x$value[match(draws, x$.draw)]
  }
  rel <- list()
  for (m in c("log-logistic", "log-normal", "weibull1", "weibull2")) {
    wf  <- beta_fixture_wf(m)
    fit <- get_brmsfit(wf)
    expect_identical(wf$meta$duration_unit, "hours")
    # relative target, every draw: exact
    res <- tls(wf, params = c("z", "ctmax"), t_ref = 60)
    rel[[m]] <- res$draws
    mid <- derive_tdt_curve(wf, temp_grid = temps)$draws
    M   <- matrix(NA_real_, max(mid$.draw), length(temps))
    M[cbind(mid$.draw, match(mid$temp, temps))] <- log10(mid$duration_model)
    want <- tls_invert_logLT(M, 0, temps)
    got  <- pick(res$draws, "CTmax", seq_len(nrow(M)))
    ok   <- is.finite(want) & is.finite(got)
    expect_gt(mean(ok), 0.9)
    expect_equal(got[ok], want[ok], tolerance = 1e-8, info = m)
    expect_equal(pick(res$draws, "z", seq_len(nrow(M))),
                 -(temps[2] - temps[1]) / (M[, 2] - M[, 1]), tolerance = 1e-8, info = m)
    # absolute target 0.3, on the first draws
    g  <- new_tdt_grid(wf, temps = t_abs, durations = dgrid)
    pp <- brms::posterior_linpred(fit, newdata = g, transform = TRUE,
                                  re_formula = NA, draw_ids = ids)
    A  <- vapply(t_abs, function(t) {
      log10(threshold_x_by_draw(pp[, g$temp == t, drop = FALSE], dgrid, 0.3))
    }, numeric(length(ids)))
    on_ct <- match(t_ct, t_abs)
    for (tr in c(60, 30)) {
      res  <- tls(wf, params = c("z", "ctmax"), t_ref = tr, target_surv = 0.3)
      want <- tls_invert_logLT(A[, on_ct], log10(tr / 60), t_ct)
      got  <- pick(res$draws, "CTmax", ids)
      ok   <- is.finite(want) & is.finite(got)
      expect_gt(mean(ok), 0.9)
      expect_lt(max(abs(got - want)[ok]), 0.02,
                label = paste(m, "t_ref", tr, "absolute CTmax"))
      if (tr == 60) za <- pick(res$draws, "z", ids)
    }
    # a draw counts only if its threshold time is on the grid at all 22
    # temperatures (tls() averages all 11 local slopes)
    slope <- (A[, match(zg + d, t_abs)] - A[, match(zg - d, t_abs)]) / (2 * d)
    zref  <- rowMeans(-1 / slope)
    ok    <- is.finite(zref) & is.finite(za)
    expect_gt(mean(ok), 0.9)
    # 0.01-0.03 on these draws; a z read off the midpoint line instead is up
    # to 0.9-1.4 away on some draws
    expect_lt(max(abs(za - zref)[ok]), 0.1, label = paste(m, "absolute z"))
  }
  # A guard that the relative path stays model-free: the swapped fixtures share
  # their draws, and only the absolute target consults the model, so relative
  # z and CTmax must be identical across models.
  for (m in names(rel)[-1]) expect_equal(rel[[m]], rel[[1]], tolerance = 1e-12, info = m)
})

# ---- gated: simulate from each model and recover it --------------------------------

# simulate_tdt()'s default design under `model` (helper-simulate.R), ready to
# fit; attr "truth" holds z and the relative-threshold CTmax at 1 h.
m_simulate <- function(model, seed) {
  raw <- simulate_tdt(seed = seed, model = model)
  std <- standardize_data(raw, temp = "T", duration = "t_hours",
                          n_total = "n_trials", n_surv = "y_alive",
                          duration_unit = "hours")
  tr  <- truth_summary(attr(raw, "truth"), model = model)
  attr(std, "truth") <- c(z = tr$z, ctmax_1h = tr$CTmax_rel_1hr)
  std
}
m_fit <- function(fun, std, ...) {
  fun(std, ..., chains = 2, iter = 800, warmup = 400, cores = 2, seed = 1,
      control = list(adapt_delta = 0.95, max_treedepth = 12), refresh = 0)
}

test_that("each new model recovers z and CTmax from data simulated under it", {
  skip_unless_brms()
  # The cached Phase 3 fits (helper-skip.R: half-decade design, seed 7); the
  # recovery tests (test-model-recovery.R) read the same fits.
  cases <- list(list(model = "log-normal"), list(model = "weibull1"),
                list(model = "weibull2"))
  for (cs in cases) {
    wf  <- load_fixture_workflow_model(cs$model)
    expect_equal(wf$meta$model, cs$model)
    tr    <- attr(wf, "truth")
    truth <- c(z = tr$z, ctmax_1h = tr$CTmax_rel_1hr)
    res <- tls(wf, params = c("z", "ctmax"), t_ref = 60)$summary
    z_hat  <- res$median[res$quantity == "z"]
    ct_hat <- res$median[res$quantity == "CTmax"]
    expect_length(z_hat, 1)
    expect_length(ct_hat, 1)
    expect_lt(abs(z_hat - truth[["z"]]) / truth[["z"]], 0.25)
    expect_lt(abs(ct_hat - truth[["ctmax_1h"]]), 1)
    ep <- brms::posterior_epred(get_brmsfit(wf), ndraws = 10)   # R-side curve evaluation
    expect_true(all(is.finite(ep)))
    pc <- predict_survival_curves(wf, temps = c(30, 33, 36))
    expect_true(is.list(pc) || is.data.frame(pc))
    # Absolute thresholds use each model's own correction. The engine's closed
    # form is checked per draw against derive_tdt_curve(), which inverts brms's
    # predicted curve numerically and knows nothing about the corrections.
    fit   <- get_brmsfit(wf)
    tbar  <- wf$meta$temp_mean
    temps <- c(30, 33, 36)
    nd    <- tls_build_grid(fit$data, by = NULL, temp = "temp_c",
                            temp_grid = temps - tbar)
    for (p in c(0.3, 0.7)) {
      eng <- suppressWarnings(tls_eval_subpars(fit, nd, wf$meta$bounds,
                                               mode = "absolute", p = p))$logLT
      num <- suppressWarnings(derive_tdt_curve(wf, temp_grid = temps,
                                               target_surv = p))$draws
      tt  <- sort(unique(num$temp))
      for (j in seq_along(temps)) {
        nj <- num[num$temp == tt[j], ]
        dd <- abs(log10(nj$duration_model) - eng[nj$.draw, j])
        dd <- dd[is.finite(dd)]
        expect_gt(length(dd), 0.9 * nrow(eng))
        expect_lt(stats::median(dd), 0.002)
        expect_lt(max(dd), 0.02)
      }
      # every reader works, and a bare brmsfit is read as the same model
      za <- tls(wf, params = "z", target_surv = p)$summary$median
      zb <- tls(fit, params = "z", target_surv = p, temp_mean = tbar)$summary$median
      expect_true(is.finite(za))
      expect_equal(zb, za, tolerance = 1e-8)
      expect_true(all(is.finite(derive_z(wf, target_surv = p)$summary$z_median)))
      # exposure in model units (these fits use hours): the temperature giving
      # a 1 h LT, about 31.5 degC; edge draws may warn (asymptotes extrapolated)
      tf <- suppressWarnings(derive_temperature_for_duration(
        wf, 1, temp_grid = seq(29, 37, 0.05), target_surv = p))
      expect_true(all(is.finite(tf$summary$temp_median)))
    }

    # Heat injury on a constant 33 degC trace (= T_bar): per draw, predicted
    # survival is the fitted curve itself (relative and absolute dose units
    # differ only by the model's offset, which cancels in survival).
    trace <- data.frame(time = seq(0, 2, by = 0.05), temp = 33)
    curve_nd <- new_tdt_grid(wf, temps = 33, durations = trace$time[-1])
    curve_pp <- posterior_linpred_tdt(wf, curve_nd)
    sp33 <- tls_eval_subpars(fit, tls_build_grid(fit$data, by = NULL,
                                                 temp = "temp_c", temp_grid = 0),
                             wf$meta$bounds, mode = "absolute", p = 0.3)
    off33 <- unname(10^(sp33$logLT - sp33$mid))
    last  <- function(x) x[x$time == max(x$time), ]
    for (shp in c("constant", "varying")) {
      rel <- predict_heat_injury(trace, wf, ndraws = 40, seed = 1, shape = shp,
                                 save_draws = TRUE)$draws
      ab  <- predict_heat_injury(trace, wf, target_surv = 0.3, ndraws = 40,
                                 seed = 1, shape = shp, save_draws = TRUE)$draws
      for (d in unique(rel$.draw)) {
        r <- rel[rel$.draw == d, ][-1, ]; a <- ab[ab$.draw == d, ][-1, ]
        expect_equal(r$survival, unname(curve_pp[d, ]), tolerance = 1e-8)
        expect_equal(a$survival, r$survival, tolerance = 1e-8)
      }
      # dose in absolute units = relative dose / 10^offset at 33 degC, where
      # the offset is the engine's absolute LT minus the midpoint, per draw
      expect_equal(last(rel)$dose / last(ab)$dose, off33[last(rel)$.draw],
                   tolerance = 1e-8, info = shp)
    }

    # planted dose: a 30 -> 36 degC ramp, truth inside the posterior band
    ramp    <- data.frame(time = seq(0, 6, by = 0.1),
                          temp = seq(30, 36, length.out = 61))
    planted <- planted_dose_from_trace(ramp, z = truth[["z"]],
                                       CTmax_1hr = truth[["ctmax_1h"]], T_c = 20)
    for (shp in c("constant", "varying")) {
      hi <- predict_heat_injury(ramp, wf, T_c = 20, ndraws = 200, seed = 1,
                                shape = shp)$summary
      end <- hi[hi$time == max(hi$time), ]
      truth_hi <- utils::tail(planted$hi_cumulative, 1)
      expect_true(end$hi_lower <= truth_hi && truth_hi <= end$hi_upper,
                  info = paste(cs$model, shp))
    }
    # relative summaries and the parameter table work, labelled with the model
    zrel <- derive_z(wf)$summary
    expect_equal(zrel$z_median, z_hat, tolerance = 1e-6)    # same z as tls()
    tab <- tdt_parameter_table(wf)
    expect_equal(attr(tab, "model"), cs$model)
    expect_true(all(is.finite(tab$median)))
  }
  # a cached file holding one model is refused when another is requested
  cache <- testthat::test_path("fixtures", "sim_fit_model_lognormal")
  expect_error(m_fit(fit_weibull, m_simulate("weibull1", seed = 7),
                     file = cache, file_refit = "never"),
               "is a log-normal model, but a weibull1 model was requested")
})

test_that("a direct log-normal fit with an absolute threshold uses inv_Phi on both sides", {
  skip_unless_brms()
  std <- m_simulate("log-normal", seed = 8)
  wf  <- m_fit(fit_lognormal, std, ctmax = ~ 1, z = ~ 1, threshold = "absolute")
  expect_equal(wf$meta$threshold, "absolute")
  lp <- brms::posterior_linpred(get_brmsfit(wf), nlpar = "mid", ndraws = 5)
  expect_true(all(is.finite(lp)))
  res <- tls(wf, params = c("z", "ctmax"), t_ref = 60)$summary
  expect_true(all(is.finite(res$median)))
  # the fitted direct line is the absolute LT50: the Stan-side inv_Phi (in
  # mid) and the R-side qnorm (in the engine) must cancel exactly
  fit <- get_brmsfit(wf)
  nd  <- tls_build_grid(fit$data, temp_grid = c(-3, 0, 3))
  eng <- tls_eval_subpars(fit, nd, wf$meta$bounds, mode = "absolute",
                          p = 0.5)$logLT
  ct  <- brms::posterior_linpred(fit, newdata = nd, nlpar = "CTmaxdev",
                                 re_formula = NA)
  lz  <- brms::posterior_linpred(fit, newdata = nd, nlpar = "logz",
                                 re_formula = NA)
  tc  <- matrix(nd$temp_c, nrow(eng), ncol(eng), byrow = TRUE)
  expect_equal(eng, wf$meta$log10_tref - (tc - ct) / exp(lz),
               tolerance = 1e-10, ignore_attr = TRUE)
  # brms computes the beta-binomial log-likelihood with extraDistr (a brms
  # Suggests); LOO is checked when it is installed.
  skip_if_not_installed("extraDistr")
  lo <- suppressWarnings(brms::loo(get_brmsfit(wf)))
  expect_s3_class(lo, "loo")
})

test_that("each model's logG is log(G) without underflow", {
  for (m in names(tdt_models())) {
    sp <- tdt_models()[[m]]
    u  <- c(-5, -1, 0, 1, 3, 6)
    expect_equal(sp$logG(u), log(sp$G(u)), tolerance = 1e-10, info = m)
    expect_equal(sp$logG(-Inf), 0, info = m)
    # far into the tail G underflows to 0, logG stays finite
    expect_true(all(is.finite(sp$logG(c(50, 200)))), info = m)
    expect_true(all(diff(sp$logG(c(-Inf, -3, 0, 3, 50, 200))) < 0), info = m)
  }
})

test_that("each model's logGc, log-scale inverses and abs_corr_dp agree with G and abs_corr", {
  for (m in names(tdt_models())) {
    sp <- tdt_models()[[m]]
    u  <- c(-2, -1, 0, 1, 2)
    expect_equal(sp$logGc(u), log1p(-sp$G(u)), tolerance = 1e-8, info = m)
    expect_equal(sp$G(0), 0.5, info = m)              # every curve's midpoint is u = 0
    # Each inverse on the side where its share is at most 1/2 (where the
    # cumulative-exposure mapping uses it), out to tails where the share
    # itself underflows to 0 on the linear scale.
    ue <- c(0, 0.5, 2, 6, 20, 50, 300)
    expect_equal(sp$Ginv_log(sp$logG(ue)), ue, tolerance = 1e-12, info = m)
    expect_equal(sp$Gcinv_log(sp$logGc(-ue)), -ue, tolerance = 1e-12, info = m)
    expect_true(all(is.finite(sp$logG(ue))) && all(is.finite(sp$logGc(-ue))), info = m)
    expect_true(all(diff(sp$logGc(c(-300, -50, -3, 0, 3))) > 0), info = m)
    # abs_corr_dp against a central difference of abs_corr in p
    low <- 0.05; up <- 0.9; k <- 6; h <- 1e-6
    for (p in c(0.1, 0.3, 0.475, 0.6, 0.85)) {
      fd <- (sp$abs_corr(low, up, k, p + h) - sp$abs_corr(low, up, k, p - h)) / (2 * h)
      expect_equal(sp$abs_corr_dp(low, up, k, p), fd, tolerance = 1e-6,
                   info = paste(m, p))
    }
  }
})
