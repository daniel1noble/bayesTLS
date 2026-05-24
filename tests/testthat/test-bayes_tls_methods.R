# Tests for print, summary, plot S3 methods on bayes_tls objects. The
# print-on-spec test runs without brms; summary and plot are gated behind
# RUN_BRMS_TESTS because they need a fitted brmsfit underneath.

test_that("print.bayes_tls handles a spec-only workflow (no fit)", {
  raw <- data.frame(
    temperature_C = rep(c(30, 32, 34), each = 3),
    exposure_h    = rep(c(1, 2, 4), times = 3),
    n             = 30L,
    alive         = c(29, 25, 5, 30, 18, 2, 28, 10, 1)
  )
  std <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                          n_total = "n", n_surv = "alive")
  wf <- fit_4pl(std, fit = FALSE)

  out <- capture.output(print(wf))
  expect_true(any(grepl("^<bayes_tls>", out)))
  expect_true(any(grepl("spec only", out)))
})

test_that("summary.bayes_tls errors when the workflow is not fitted", {
  raw <- data.frame(
    temperature_C = rep(c(30, 32, 34), each = 3),
    exposure_h    = rep(c(1, 2, 4), times = 3),
    n             = 30L,
    alive         = c(29, 25, 5, 30, 18, 2, 28, 10, 1)
  )
  std <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                          n_total = "n", n_surv = "alive")
  wf <- fit_4pl(std, fit = FALSE)
  expect_error(summary(wf), "no fit")
})

test_that("plot.bayes_tls errors when the workflow is not fitted", {
  raw <- data.frame(
    temperature_C = rep(c(30, 32, 34), each = 3),
    exposure_h    = rep(c(1, 2, 4), times = 3),
    n             = 30L,
    alive         = c(29, 25, 5, 30, 18, 2, 28, 10, 1)
  )
  std <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                          n_total = "n", n_surv = "alive")
  wf <- fit_4pl(std, fit = FALSE)
  expect_error(plot(wf), "no fit")
})

test_that("get_brmsfit errors on a spec-only (unfitted) workflow", {
  raw <- data.frame(
    temperature_C = rep(c(30, 32, 34), each = 3),
    exposure_h    = rep(c(1, 2, 4), times = 3),
    n             = 30L,
    alive         = c(29, 25, 5, 30, 18, 2, 28, 10, 1)
  )
  std <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                          n_total = "n", n_surv = "alive")
  wf <- fit_4pl(std, fit = FALSE)
  expect_false(has_fit(wf))
  expect_error(get_brmsfit(wf), "no fit")
  # Also rejects a non-workflow object outright.
  expect_error(get_brmsfit(list(fit = NULL)), "no fit")
})

test_that("get_brmsfit returns the underlying brmsfit on a fitted workflow", {
  skip_unless_brms()

  wf  <- load_fixture_workflow()
  fit <- get_brmsfit(wf)
  expect_s3_class(fit, "brmsfit")
  expect_identical(fit, wf$fit)
  # The returned object is usable by brms helpers downstream.
  expect_true(brms::ndraws(fit) > 0)
})

test_that("summary.bayes_tls delegates to brms and forwards `...`", {
  skip_unless_brms()

  wf <- load_fixture_workflow()
  s  <- summary(wf)

  # The method delegates to brms::summary.brmsfit, so the object is the brms
  # summary, not a custom structure.
  expect_s3_class(s, "brmssummary")
  # The population-level coefficient table carries our 4PL sub-parameters
  # (brms strips the `b_` prefix in the printed/fixed table).
  expect_true(all(c("Estimate", "Rhat", "Bulk_ESS") %in% colnames(s$fixed)))
  expect_true("mid_Intercept" %in% rownames(s$fixed))

  # The fixture is a small, well-behaved fit: Rhat under a safe ceiling.
  expect_lt(max(s$fixed[, "Rhat"], na.rm = TRUE), 1.1)

  # `...` is forwarded to brms::summary.brmsfit — e.g. credible-interval width.
  s90 <- summary(wf, prob = 0.9)
  expect_true(any(grepl("90% CI", colnames(s90$fixed))))
})

test_that("printing summary(bayes_tls) shows the brms summary of our model", {
  skip_unless_brms()

  wf  <- load_fixture_workflow()
  out <- capture.output(print(summary(wf)))
  # brms's summary print header plus our model's parameter names.
  expect_true(any(grepl("Family:", out)))
  expect_true(any(grepl("mid_Intercept", out)))
})

test_that("plot.bayes_tls delegates to brms and returns its trace-plot output", {
  skip_unless_brms()

  wf <- load_fixture_workflow()
  g  <- plot(wf)
  # brms::plot.brmsfit returns the bayesplot grid object(s) it draws.
  expect_type(g, "list")
  expect_s3_class(g[[1]], "bayesplot_grid")

  # `...` is forwarded: subset to a single parameter via brms's `variable`.
  g2 <- plot(wf, variable = "b_mid_Intercept")
  expect_type(g2, "list")

  # A non-existent parameter raises brms's own error.
  expect_error(plot(wf, variable = "nonexistent_param"),
               "missing in the draws object")
})
