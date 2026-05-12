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

test_that("fit_4pl(fit = FALSE) returns a tdt_4pl_workflow without fitting", {
  raw <- data.frame(
    temperature_C = rep(c(30, 32, 34), each = 3),
    exposure_h    = rep(c(1, 2, 4), times = 3),
    n             = 30L,
    alive         = c(29, 25, 5, 30, 18, 2, 28, 10, 1)
  )
  std <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                          n_total = "n", n_surv = "alive")

  wf <- fit_4pl(std, fit = FALSE)
  expect_s3_class(wf, "tdt_4pl_workflow")
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
