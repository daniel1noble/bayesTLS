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

test_that("standardize_data respects an explicit temp_mean", {
  raw <- make_raw()
  std <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                          n_total = "n", n_surv = "alive",
                          temp_mean = 31)
  expect_equal(attr(std, "tdt_meta")$temp_mean, 31)
  expect_equal(std$temp_c, std$temp - 31)
})
