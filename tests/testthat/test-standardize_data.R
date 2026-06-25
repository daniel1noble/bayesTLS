make_raw <- function() {
  data.frame(
    temperature_C = rep(c(30, 32, 34), each = 4),
    exposure_h    = rep(c(1, 2, 4, 8), times = 3),
    n             = 30L,
    alive         = c(29, 28, 25, 5, 30, 27, 18, 2, 28, 22, 10, 1),
    Tank          = rep(c("T1", "T2", "T3"), times = 4),
    Date          = rep(c("D1", "D2"),        times = 6)
  )
}

test_that("standardize_data produces the standard column set with n_surv", {
  raw <- make_raw()
  std <- standardize_data(raw,
                          temp = "temperature_C", duration = "exposure_h",
                          n_total = "n", n_surv = "alive")
  expect_s3_class(std, "tbl_df")
  expect_true(all(c("temp", "duration", "logd", "temp_c",
                    "n_total", "n_surv", "n_dead", "survival") %in% names(std)))
  expect_equal(nrow(std), nrow(raw))
  expect_true(all(std$n_surv + std$n_dead == std$n_total))
  expect_equal(std$temp_c, std$temp - mean(std$temp))
})

test_that("standardize_data handles each count specification consistently", {
  raw <- make_raw()
  raw$mortality_prop <- (raw$n - raw$alive) / raw$n
  raw$n_dead         <- raw$n - raw$alive
  raw$survival_prop  <- raw$alive / raw$n

  by_surv <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                              n_total = "n", n_surv = "alive")
  by_dead <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                              n_total = "n", n_dead = "n_dead")
  by_prop <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                              n_total = "n", survival = "survival_prop")
  by_mort <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                              n_total = "n", mortality = "mortality_prop")

  expect_equal(by_surv$n_surv, by_dead$n_surv)
  expect_equal(by_surv$n_surv, by_prop$n_surv)
  expect_equal(by_surv$n_surv, by_mort$n_surv)
})

test_that("standardize_data rejects no/multiple count specs", {
  raw <- make_raw()
  expect_error(standardize_data(raw, temp = "temperature_C",
                                duration = "exposure_h", n_total = "n"),
               "exactly one of")
  expect_error(standardize_data(raw, temp = "temperature_C",
                                duration = "exposure_h",
                                n_total = "n", n_surv = "alive", n_dead = "n"),
               "exactly one of")
})

test_that("standardize_data errors on missing input columns", {
  raw <- make_raw()
  expect_error(standardize_data(raw, temp = "missing_col",
                                duration = "exposure_h",
                                n_total = "n", n_surv = "alive"),
               "Missing input columns")
})

test_that("standardize_data attaches metadata with the right fields", {
  raw <- make_raw()
  std <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                          n_total = "n", n_surv = "alive",
                          random_effects = c("Date", "Tank"),
                          duration_unit = "hours")
  meta <- attr(std, "tdt_meta")
  expect_equal(meta$duration_unit,  "hours")
  expect_equal(meta$random_effects, c("Date", "Tank"))
  expect_equal(meta$temp_mean, mean(std$temp))
  # Random-effect columns should be coerced to factors
  expect_s3_class(std$Tank, "factor")
  expect_s3_class(std$Date, "factor")
})

test_that("standardize_data proportion path produces a Beta-ready frame", {
  raw <- make_raw()
  raw$fvfm <- c(0, 0.2, 0.5, 0.95, 1, 0.3, 0.6, 0.8, 0.1, 0.4, 0.7, 0.9)
  std <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                          proportion = "fvfm")
  meta <- attr(std, "tdt_meta")
  expect_equal(meta$response_type, "proportion")
  expect_equal(meta$response_var,  "survival")
  expect_false("n_total" %in% names(std))
  expect_false("n_surv"  %in% names(std))
  expect_true(all(c("temp", "duration", "logd", "temp_c", "survival") %in%
                  names(std)))
  # Exact 0 and 1 are clamped into the open interval for the Beta likelihood.
  expect_true(all(std$survival > 0 & std$survival < 1))
  expect_equal(min(std$survival), 0.001)
  expect_equal(max(std$survival), 0.999)
})

test_that("proportion path records proportion_eps in metadata; count path does not", {
  raw <- make_raw()
  raw$fvfm <- c(0, 0.2, 0.5, 0.95, 1, 0.3, 0.6, 0.8, 0.1, 0.4, 0.7, 0.9)
  std_p <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                            proportion = "fvfm", proportion_eps = 0.005)
  expect_equal(attr(std_p, "tdt_meta")$proportion_eps, 0.005)

  std_c <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                            n_total = "n", n_surv = "alive")
  expect_null(attr(std_c, "tdt_meta")$proportion_eps)
})

test_that("proportion path warns when overwriting an unrelated `survival` column", {
  raw <- make_raw()
  raw$fvfm     <- c(0.1, 0.2, 0.5, 0.95, 0.4, 0.3, 0.6, 0.8, 0.1, 0.4, 0.7, 0.9)
  raw$survival <- "unrelated"                  # pre-existing, not the source
  expect_warning(
    std <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                            proportion = "fvfm"),
    "overwrites the existing `survival` column")
  expect_equal(std$survival, pmin(pmax(raw$fvfm, 0.001), 0.999))

  # No warning when the proportion source *is* `survival`.
  raw2 <- make_raw()
  raw2$survival <- c(0.1, 0.2, 0.5, 0.95, 0.4, 0.3, 0.6, 0.8, 0.1, 0.4, 0.7, 0.9)
  expect_no_warning(
    standardize_data(raw2, temp = "temperature_C", duration = "exposure_h",
                     proportion = "survival"))
})

test_that("standardize_data proportion path needs no n_total and rejects mixing", {
  raw <- make_raw()
  raw$fvfm <- 0.5
  # No n_total required for a proportion response.
  expect_s3_class(
    standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                     proportion = "fvfm"),
    "tbl_df"
  )
  # Proportion combined with a count spec is rejected.
  expect_error(
    standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                     n_total = "n", n_surv = "alive", proportion = "fvfm"),
    "exactly one of"
  )
})

test_that("standardize_data count path still requires n_total", {
  raw <- make_raw()
  expect_error(
    standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                     n_surv = "alive"),
    "require `n_total`"
  )
})

test_that("standardize_data respects an explicit temp_mean", {
  raw <- make_raw()
  std <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                          n_total = "n", n_surv = "alive",
                          temp_mean = 31)
  expect_equal(attr(std, "tdt_meta")$temp_mean, 31)
  expect_equal(std$temp_c, std$temp - 31)
})

test_that("counts passed to survival=/mortality= error (not silently clamped)", {
  raw <- data.frame(t = c(30, 32, 34), d = c(1, 1, 1), n = 30,
                    s = c(29, 18, 2))            # survivor COUNTS, not proportions
  expect_error(
    standardize_data(raw, temp = "t", duration = "d", n_total = "n", survival = "s"),
    "values > 1")
  rawm <- data.frame(t = c(30, 32, 34), d = 1, n = 30, m = c(1, 12, 28))
  expect_error(
    standardize_data(rawm, temp = "t", duration = "d", n_total = "n", mortality = "m"),
    "values > 1")
})

test_that("survivor counts above n_total warn (clamped, not silently fabricated)", {
  raw <- data.frame(t = c(30, 32), d = c(1, 1), n = c(10, 10),
                    alive = c(12, 8))               # 12 > 10 is a data-entry error
  expect_warning(
    d <- standardize_data(raw, temp = "t", duration = "d",
                          n_total = "n", n_surv = "alive"),
    "outside \\[0, n_total\\]")
  expect_true(all(d$n_surv <= d$n_total))           # clamped to a valid range
})

test_that("a continuous proportion outside [0, 1] errors", {
  raw <- data.frame(t = c(30, 32, 34), d = 1, p = c(0.9, 0.5, 1.4))
  expect_error(
    standardize_data(raw, temp = "t", duration = "d", proportion = "p"),
    "outside \\[0, 1\\]")
})

test_that("a single assay temperature warns about an unidentified z", {
  raw <- data.frame(t = 30, d = c(1, 2, 4), n = 30, alive = c(29, 20, 5))
  expect_warning(
    standardize_data(raw, temp = "t", duration = "d", n_total = "n", n_surv = "alive"),
    "not identified")
})

test_that("clobbering a transform-output column (temp_c) warns", {
  raw <- data.frame(temperature_C = c(30, 32, 34), exposure_h = 1, n = 30,
                    alive = c(29, 20, 5), temp_c = c("a", "b", "c"))  # unrelated column
  expect_warning(
    standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                     n_total = "n", n_surv = "alive"),
    "overwrites existing column")
})
