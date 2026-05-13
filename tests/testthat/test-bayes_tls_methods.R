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

test_that("summary.bayes_tls returns the expected structure", {
  skip_unless_brms()

  wf <- load_fixture_workflow()
  s  <- summary(wf)

  expect_s3_class(s, "summary.bayes_tls")
  expect_named(s, c("meta", "shape", "parameters", "diagnostics"),
               ignore.order = TRUE)
  expect_true(all(c("variable", "mean", "rhat", "ess_bulk") %in%
                  names(s$parameters)))
  expect_named(s$diagnostics,
               c("max_rhat", "min_ess_bulk", "min_ess_tail",
                 "divergences", "treedepth_hits"),
               ignore.order = TRUE)

  # The fixture is a small, well-behaved fit: expect Rhat under a safe ceiling.
  expect_lt(s$diagnostics$max_rhat, 1.1)
})

test_that("print.summary.bayes_tls prints diagnostics + posterior table", {
  skip_unless_brms()

  wf  <- load_fixture_workflow()
  out <- capture.output(print(summary(wf)))
  expect_true(any(grepl("^<bayes_tls summary>", out)))
  expect_true(any(grepl("HMC diagnostics", out)))
  expect_true(any(grepl("Posterior summary", out)))
  expect_true(any(grepl("b_mid_Intercept", out)))
})

test_that("plot.bayes_tls returns a ggplot of trace plots", {
  skip_unless_brms()

  wf <- load_fixture_workflow()
  g  <- plot(wf)
  expect_s3_class(g, "ggplot")

  # Subset to a single parameter — should still return a ggplot.
  g2 <- plot(wf, pars = "b_mid_Intercept")
  expect_s3_class(g2, "ggplot")

  expect_error(plot(wf, pars = "nonexistent_param"), "Unknown parameter")
})
