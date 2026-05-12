small_std <- function() {
  raw <- data.frame(
    temperature_C = rep(c(30, 32, 34), each = 4),
    exposure_h    = rep(c(1, 2, 4, 8), times = 3),
    n             = 30L,
    alive         = c(29, 28, 25, 5, 30, 27, 18, 2, 28, 22, 10, 1)
  )
  standardize_data(raw,
                   temp = "temperature_C", duration = "exposure_h",
                   n_total = "n", n_surv = "alive")
}

test_that("make_4pl_priors returns an Intercept + general b prior for each sub-param", {
  p <- make_4pl_priors(small_std())
  expect_s3_class(p, "brmsprior")
  nlpars <- p$nlpar
  for (nm in c("lowraw", "upraw", "logk", "mid")) {
    rows <- p[p$nlpar == nm & p$class == "b", ]
    expect_true(any(rows$coef == "Intercept"),
                info = paste("Intercept prior for", nm))
    expect_true(any(rows$coef == ""),
                info = paste("general b prior for", nm))
  }
  expect_true(any(p$class == "phi"))
})

test_that("make_4pl_priors random_effects adds one sd prior per group", {
  p   <- make_4pl_priors(small_std(),
                         random_effects = c("Date", "Tank"))
  sds <- p[p$class == "sd", ]
  expect_equal(nrow(sds), 2)
  expect_setequal(sds$group, c("Date", "Tank"))
})

test_that("make_4pl_priors shifts Intercept centres for non-default bounds (PSII)", {
  std <- small_std()
  p_default <- make_4pl_priors(std)
  p_psii    <- make_4pl_priors(std, lower = 0.85, upper = 1)

  parse_mean <- function(prior_str) {
    as.numeric(sub("normal\\(([^,]+),.*", "\\1", prior_str))
  }
  m_low_default <- parse_mean(
    p_default[p_default$nlpar == "lowraw" & p_default$coef == "Intercept", "prior"]
  )
  m_low_psii <- parse_mean(
    p_psii[p_psii$nlpar == "lowraw" & p_psii$coef == "Intercept", "prior"]
  )
  # PSII bounds shift lowraw centre more negative (low asymptote sits closer
  # to the lower edge of the (0.85, 0.92) interval than 0.02 sits in
  # (0.001, 0.499)).
  expect_lt(m_low_psii, m_low_default)
})
