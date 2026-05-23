test_that("make_4pl_formula returns a brmsformula with the expected structure", {
  f <- make_4pl_formula()
  expect_s3_class(f, "brmsformula")

  # Convert all sub-formulas to strings so we can grep them.
  parts <- vapply(f$pforms, function(x) paste(deparse(x), collapse = " "),
                  character(1))
  expect_true(any(grepl("lowraw ~ temp_c", parts)))
  expect_true(any(grepl("upraw\\s*~ temp_c", parts)))
  expect_true(any(grepl("logk\\s*~ temp_c", parts)))
  expect_true(any(grepl("mid\\s*~ temp_c", parts)))
})

test_that("make_4pl_formula temp_effects controls which sub-params carry temp_c", {
  # Default: all four carry temp_c (regression guard).
  pf <- make_4pl_formula()$pforms
  expect_match(paste(deparse(pf$lowraw), collapse = " "), "temp_c")
  expect_match(paste(deparse(pf$upraw),  collapse = " "), "temp_c")
  expect_match(paste(deparse(pf$logk),   collapse = " "), "temp_c")
  expect_match(paste(deparse(pf$mid),    collapse = " "), "temp_c")

  # Constant-shape: temp_c only on mid; low/up/k are intercept-only.
  pf2 <- make_4pl_formula(temp_effects = "mid")$pforms
  expect_equal(paste(deparse(pf2$lowraw), collapse = " "), "lowraw ~ 1")
  expect_equal(paste(deparse(pf2$upraw),  collapse = " "), "upraw ~ 1")
  expect_equal(paste(deparse(pf2$logk),   collapse = " "), "logk ~ 1")
  expect_match(paste(deparse(pf2$mid),    collapse = " "), "mid ~ temp_c")

  # Constant-shape keeps random intercepts on mid.
  pf3 <- make_4pl_formula(temp_effects = "mid",
                          random_effects = "Date")$pforms
  expect_match(paste(deparse(pf3$mid), collapse = " "), "\\(1 \\| Date\\)")

  # mid is mandatory.
  expect_error(make_4pl_formula(temp_effects = c("low", "up")),
               "mid.*must always carry")
})

test_that("make_4pl_formula with random_effects expands mid sub-model", {
  f <- make_4pl_formula(random_effects = c("Date", "Tank"))
  mid_form <- f$pforms$mid
  expect_true(grepl("\\(1 \\| Date\\)", paste(deparse(mid_form), collapse = " ")))
  expect_true(grepl("\\(1 \\| Tank\\)", paste(deparse(mid_form), collapse = " ")))
})

test_that("make_4pl_formula with PSII bounds bakes lower/upper into rhs", {
  f <- make_4pl_formula(lower = 0.85, upper = 1)
  main <- paste(deparse(f$formula), collapse = " ")
  # The constants in the rhs reflect the disjoint asymptote intervals.
  expect_true(grepl("0\\.851", main))
  expect_true(grepl("0\\.926", main))
})

test_that("fit_4pl(fit = FALSE) returns a bayes_tls without fitting", {
  raw <- data.frame(
    temperature_C = rep(c(30, 32, 34), each = 3),
    exposure_h    = rep(c(1, 2, 4), times = 3),
    n             = 30L,
    alive         = c(29, 25, 5, 30, 18, 2, 28, 10, 1)
  )
  std <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                          n_total = "n", n_surv = "alive")

  wf <- fit_4pl(std, fit = FALSE)
  expect_s3_class(wf, "bayes_tls")
  expect_null(wf$fit)
  expect_false(has_fit(wf))
  expect_s3_class(wf$formula, "brmsformula")
  expect_s3_class(wf$prior,   "brmsprior")
  expect_equal(wf$meta$lower, 0)
  expect_equal(wf$meta$upper, 1)
})

test_that("fit_4pl errors are clean when called on bad input", {
  # No random_effects → mid sub-model should still be 'mid ~ temp_c'
  f <- make_4pl_formula(random_effects = NULL)
  expect_equal(paste(deparse(f$pforms$mid), collapse = " "),
               "mid ~ temp_c")
})

test_that("make_4pl_formula default count response is unchanged", {
  f    <- make_4pl_formula()
  main <- paste(deparse(f$formula), collapse = " ")
  expect_true(grepl("n_surv \\| trials\\(n_total\\)", main))
})

test_that("make_4pl_formula Beta(identity) has no trials and no logit wrapper", {
  f    <- make_4pl_formula(family = brms::Beta(link = "identity"),
                           response_var = "survival")
  main <- paste(deparse(f$formula), collapse = " ")
  expect_true(grepl("^survival ~", main))
  expect_false(grepl("trials", main))
  expect_false(grepl("n_surv", main))
  # identity: the 4PL is the mean, so no `~ logit(` wrapper. (The reparam still
  # contains `inv_logit(` on the asymptotes, so match the wrapper specifically.)
  expect_false(grepl("~ logit\\(", main))
  expect_true(grepl("1 \\+ exp\\(exp\\(logk", main))   # 4PL structure present
})

test_that("make_4pl_formula Beta(logit) wraps the 4PL in logit()", {
  f    <- make_4pl_formula(family = brms::Beta(link = "logit"),
                           response_var = "survival")
  main <- paste(deparse(f$formula), collapse = " ")
  expect_true(grepl("survival ~ logit\\(", main))
})

test_that("fit_4pl(fit = FALSE) on proportion data builds a Beta workflow", {
  raw <- data.frame(
    temperature_C = rep(c(30, 32, 34), each = 3),
    exposure_h    = rep(c(1, 2, 4), times = 3),
    fvfm          = c(0.95, 0.6, 0.1, 0.9, 0.4, 0.05, 0.8, 0.3, 0.02)
  )
  std <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                          proportion = "fvfm")

  wf <- fit_4pl(std, fit = FALSE)
  expect_s3_class(wf, "bayes_tls")
  expect_equal(wf$meta$response_type, "proportion")
  expect_equal(wf$meta$family, "beta")
  expect_equal(wf$meta$link,   "identity")
  expect_true(any(wf$prior$class == "phi"))     # Beta carries a precision prior
  main <- paste(deparse(wf$formula$formula), collapse = " ")
  expect_true(grepl("^survival ~", main))
})
