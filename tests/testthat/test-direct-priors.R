# Tier-1 fast, deterministic, brms-free tests for make_4pl_priors() direct mode
# and the fit_4pl() direct wiring. No model is fitted.

# A small standardised-shape data frame with a grouping factor.
dd <- data.frame(
  logd       = log10(rep(c(1, 2, 4, 8), 3)),
  temp_c     = rep(c(-2, 0, 2), each = 4),
  n_surv     = 5L, n_total = 10L,
  life_stage = factor(rep(c("a", "b", "c"), each = 4))
)
ps <- function(p, nlpar, coef = "", class = "b") {
  d <- as.data.frame(p)
  v <- d$prior[d$nlpar == nlpar & d$coef == coef & d$class == class]
  if (length(v)) v[1] else NA_character_
}

test_that("midpoint priors are byte-for-byte unchanged (frozen literals)", {
  p <- make_4pl_priors(dd)
  expect_equal(ps(p, "lowraw", "Intercept"), "normal(-3.227262, 1)")
  expect_equal(ps(p, "upraw",  "Intercept"), "normal(3.227262, 1)")
  expect_equal(ps(p, "logk",   "Intercept"), "normal(0.693147, 1)")
  expect_equal(ps(p, "lowraw"), "normal(0, 0.5)")   # general slope prior
  expect_equal(ps(p, "logk"),   "normal(0, 0.3)")
  expect_true(any(as.data.frame(p)$class == "phi"))
  expect_false(any(c("CTmaxdev", "logz") %in% as.data.frame(p)$nlpar))
})

test_that("supplying ctmax/z switches to direct priors (CTmaxdev + logz, no mid)", {
  p <- make_4pl_priors(dd, ctmax = ~ 1, z = ~ 1)
  nl <- as.data.frame(p)$nlpar
  expect_true(all(c("CTmaxdev", "logz") %in% nl))
  expect_false("mid" %in% nl)
  # The level prior lands on the Intercept; a separate mean-zero prior covers
  # contrasts/slopes (so per-group z/CTmax are coding-invariant).
  expect_equal(ps(p, "CTmaxdev", "Intercept"), "normal(0.000000, 10)")
  expect_equal(ps(p, "CTmaxdev"),              "normal(0, 10)")
  expect_equal(ps(p, "logz", "Intercept"),     "normal(1.098612, 0.7)")
  expect_equal(ps(p, "logz"),                  "normal(0, 0.7)")
})

test_that("CTmaxdev/logz contrast prior is mean-zero -> coding-invariant", {
  # Treatment coding (~ G): Intercept carries the level centre; the between-group
  # contrast must be centred on ZERO, not on log(3) (the bug that biased the
  # between-group z ratio toward 3x).
  p_trt <- make_4pl_priors(dd, ctmax = ~ life_stage, z = ~ life_stage)
  expect_equal(ps(p_trt, "logz", "Intercept"),     "normal(1.098612, 0.7)")
  expect_equal(ps(p_trt, "logz"),                  "normal(0, 0.7)")  # -> contrasts
  expect_equal(ps(p_trt, "CTmaxdev", "Intercept"), "normal(0.000000, 10)")
  expect_equal(ps(p_trt, "CTmaxdev"),              "normal(0, 10)")

  # Cell-means coding (~ 0 + G): every level carries the same centre, so there is
  # no reference group to bias.
  d  <- as.data.frame(make_4pl_priors(dd, ctmax = ~ 0 + life_stage,
                                      z = ~ 0 + life_stage))
  lz <- d$prior[d$nlpar == "logz" & d$coef != ""]
  expect_true(length(lz) == 3L && all(lz == "normal(1.098612, 0.7)"))
})

test_that("centred asymptote priors are byte-identical between midpoint and direct", {
  pm <- make_4pl_priors(dd)
  pd <- make_4pl_priors(dd, ctmax = ~ 1, z = ~ 1)
  for (nlp in c("lowraw", "upraw", "logk")) {
    expect_identical(ps(pd, nlp, "Intercept"), ps(pm, nlp, "Intercept"))
    expect_identical(ps(pd, nlp), ps(pm, nlp))   # general slope prior too
  }
})

test_that("cell-means ctmax/z yields one centred asymptote prior per factor level", {
  p <- make_4pl_priors(dd, ctmax = ~ 0 + life_stage, z = ~ 0 + life_stage)
  d <- as.data.frame(p)
  low_coefs <- d$coef[d$nlpar == "lowraw" & d$coef != ""]
  expect_setequal(low_coefs, c("life_stagea", "life_stageb", "life_stagec"))
  # every per-level asymptote prior carries the centred value, not the flat slope
  expect_true(all(d$prior[d$nlpar == "lowraw" & d$coef != ""] == "normal(-3.227262, 1)"))
})

test_that("random effects attach to CTmaxdev/logz only, never to the shape", {
  p  <- make_4pl_priors(dd, ctmax = ~ 0 + life_stage + (1 | batch),
                        z = ~ 0 + life_stage + (1 | batch))
  sd <- as.data.frame(p)[as.data.frame(p)$class == "sd", ]
  expect_setequal(sd$nlpar, c("CTmaxdev", "logz"))
  expect_true(all(sd$group == "batch"))
  expect_false(any(sd$nlpar %in% c("lowraw", "upraw", "logk")))
})

test_that("phi prior present for beta-binomial, omitted when prior_phi = NULL", {
  expect_true(any(as.data.frame(make_4pl_priors(dd, ctmax = ~ 1))$class == "phi"))
  expect_false(any(as.data.frame(
    make_4pl_priors(dd, ctmax = ~ 1, prior_phi = NULL))$class == "phi"))
})

test_that("bare cell-means (0 + G) emits no redundant global 'b' prior", {
  # Every coefficient already has a per-level prior, so the global mean-zero
  # `b` prior would be redundant and brms warns about it (CI = error-on=warning).
  p <- as.data.frame(make_4pl_priors(dd, ctmax = ~ 0 + life_stage,
                                     z = ~ 0 + life_stage))
  for (nlp in c("CTmaxdev", "logz", "lowraw", "upraw", "logk")) {
    has_global <- any(p$nlpar == nlp & p$coef == "" & p$class == "b")
    expect_false(has_global)
  }
})

test_that("cell-means with extra terms keeps the global mean-zero 'b' prior", {
  # An interaction/slope term is not covered by the per-level priors, so the
  # global mean-zero `b` prior must remain.
  p <- as.data.frame(make_4pl_priors(
    dd, ctmax = ~ 0 + life_stage + temp_c:life_stage, z = ~ 0 + life_stage))
  expect_true(any(p$nlpar == "CTmaxdev" & p$coef == "" & p$class == "b"))
  # treatment coding still carries the global prior (covers contrasts)
  pt <- as.data.frame(make_4pl_priors(dd, ctmax = ~ life_stage, z = ~ life_stage))
  expect_true(any(pt$nlpar == "logz" & pt$coef == "" & pt$class == "b"))
})

test_that("direct prior object is complete/valid for brms (every coef covered)", {
  f <- make_4pl_formula(ctmax = ~ 0 + life_stage, z = ~ 0 + life_stage)
  p <- make_4pl_priors(dd, ctmax = ~ 0 + life_stage, z = ~ 0 + life_stage)
  expect_no_error(
    brms::validate_prior(p, f, data = dd,
                         family = brms::beta_binomial(link = "identity")))
})

test_that("fit_4pl(fit = FALSE) records the direct parameterisation in meta", {
  d <- dd; attr(d, "tdt_meta") <- list(temp_mean = 35, duration_unit = "hours",
                                       response_type = "count")
  wf <- fit_4pl(d, ctmax = ~ 0 + life_stage, z = ~ 0 + life_stage, fit = FALSE)
  expect_equal(wf$meta$parameterization, "direct")
  expect_equal(wf$meta$threshold, "relative")
  expect_true("CTmaxdev" %in% names(wf$formula$pforms))
  # midpoint stays midpoint
  wm <- fit_4pl(d, fit = FALSE)
  expect_equal(wm$meta$parameterization, "midpoint")
  expect_true("mid" %in% names(wm$formula$pforms))
})

test_that("fit_4pl normalises single-factor treatment coding to cell-means", {
  d <- dd; attr(d, "tdt_meta") <- list(temp_mean = 0, duration_unit = "minutes",
                                       response_type = "count")
  wt <- fit_4pl(d, ctmax = ~ life_stage,     z = ~ life_stage,     fit = FALSE)
  wc <- fit_4pl(d, ctmax = ~ 0 + life_stage, z = ~ 0 + life_stage, fit = FALSE)
  # treatment (~ G) and cell-means (~ 0 + G) build the identical model ...
  for (nl in c("CTmaxdev", "logz", "lowraw", "upraw", "logk"))
    expect_identical(deparse(wt$formula$pforms[[nl]]),
                     deparse(wc$formula$pforms[[nl]]))
  # ... with identical priors ...
  ord <- function(p) { x <- as.data.frame(p); x[order(x$nlpar, x$coef, x$class), ] }
  expect_equal(ord(wt$prior), ord(wc$prior), ignore_attr = TRUE)
  # ... and per-level CTmax/z priors (cell-means), not Intercept + contrasts.
  pr <- as.data.frame(wt$prior)
  expect_false(any(pr$nlpar == "logz" & pr$coef == "Intercept"))
  expect_true(all(c("life_stagea", "life_stageb", "life_stagec") %in%
                  pr$coef[pr$nlpar == "logz"]))
})

test_that("t_ref -> log10_tref respects the data's time unit", {
  d_h <- dd; attr(d_h, "tdt_meta") <- list(temp_mean = 35, duration_unit = "hours",
                                           response_type = "count")
  d_m <- dd; attr(d_m, "tdt_meta") <- list(temp_mean = 35, duration_unit = "minutes",
                                           response_type = "count")
  # 1-hour reference (t_ref = 60 min): hours model -> log10(1) = 0; minutes -> log10(60)
  expect_equal(fit_4pl(d_h, ctmax = ~ 1, fit = FALSE)$meta$log10_tref, 0)
  expect_equal(fit_4pl(d_m, ctmax = ~ 1, fit = FALSE)$meta$log10_tref, log10(60))
})
