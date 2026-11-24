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
  expect_silent(tdt_require_loglogistic(wf))
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

test_that("absolute thresholds and heat-injury prediction are refused for new models until supported", {
  std <- m_std()
  wf <- fit_lognormal(std, fit = FALSE)
  expect_error(tls(wf, target_surv = "absolute"), "not yet available for the log-normal")
  expect_error(predict_heat_injury(data.frame(time = 0:1, temp = c(30, 30)), wf),
               "log-logistic fits only")
})

# ---- gated: simulate from each model and recover it --------------------------------

m_simulate <- function(model, seed, truth = list(ell = 0.05, u = 0.95, k = 6,
                                                 m_beta0 = 1.5, m_beta1 = -0.18,
                                                 T_bar = 33)) {
  set.seed(seed)
  G <- tdt_model_spec(model)$G
  d <- expand.grid(T = c(30, 33, 36), t_hours = c(0.05, 0.5, 5, 50), rep = 1:5)
  d$log10_t_min <- log10(d$t_hours * 60)
  d$mid_t  <- truth$m_beta0 + truth$m_beta1 * (d$T - truth$T_bar)
  d$p_true <- truth$ell + (truth$u - truth$ell) * G(truth$k * (d$log10_t_min - d$mid_t))
  d$n_trials <- 30L
  d$p_draw   <- rbeta(nrow(d), d$p_true * 8, (1 - d$p_true) * 8)
  d$y_alive  <- rbinom(nrow(d), size = d$n_trials, prob = d$p_draw)
  std <- standardize_data(d, temp = "T", duration = "t_hours",
                          n_total = "n_trials", n_surv = "y_alive",
                          duration_unit = "hours")
  attr(std, "truth") <- c(z = -1 / truth$m_beta1,
                          ctmax_1h = truth$T_bar + (log10(60) - truth$m_beta0) / truth$m_beta1)
  std
}
m_fit <- function(fun, std, ...) {
  fun(std, ..., chains = 2, iter = 800, warmup = 400, cores = 2, seed = 1,
      control = list(adapt_delta = 0.95, max_treedepth = 12), refresh = 0)
}

test_that("each new model recovers z and CTmax from data simulated under it", {
  skip_unless_brms()
  cache <- tempfile("bayesTLS-models-")        # shared by the cache check below
  cases <- list(list(model = "log-normal", fun = fit_lognormal, args = list()),
                list(model = "weibull1",   fun = fit_weibull,   args = list(type = 1)),
                list(model = "weibull2",   fun = fit_weibull,   args = list(type = 2)))
  for (cs in cases) {
    std <- m_simulate(cs$model, seed = 7)
    wf  <- do.call(m_fit, c(list(cs$fun, std), cs$args,
                            if (cs$model == "log-normal") list(file = cache)))
    expect_equal(wf$meta$model, cs$model)
    truth <- attr(std, "truth")
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
    # every reader that inverts at an absolute threshold refuses, not just tls()
    msg <- paste("not yet available for the", cs$model)
    expect_error(tls(wf, target_surv = "absolute"), msg)
    expect_error(tls(get_brmsfit(wf), target_surv = 0.5), msg)  # bare brmsfit
    expect_error(derive_z(wf, target_surv = "absolute"), msg)
    expect_error(derive_temperature_for_duration(wf, 60, temp_grid = 30:36,
                                                 target_surv = 0.5), msg)
    expect_error(predict_heat_injury(data.frame(time = 0:1, temp = c(33, 33)), wf),
                 "log-logistic fits only")
    # relative summaries and the parameter table work, labelled with the model
    zrel <- derive_z(wf)$summary
    expect_equal(zrel$z_median, z_hat, tolerance = 1e-6)    # same z as tls()
    tab <- tdt_parameter_table(wf)
    expect_equal(attr(tab, "model"), cs$model)
    expect_true(all(is.finite(tab$median)))
  }
  # a cached file holding one model is refused when another is requested
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
  # brms computes the beta-binomial log-likelihood with extraDistr (a brms
  # Suggests); LOO is checked when it is installed.
  skip_if_not_installed("extraDistr")
  lo <- suppressWarnings(brms::loo(get_brmsfit(wf)))
  expect_s3_class(lo, "loo")
})
