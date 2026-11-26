# Censored continuous-proportion readings: standardize_data(detection_limit =)
# flags readings beyond an instrument limit, and the fitters add brms's cens()
# term so those readings contribute a probability instead of a density. Fast
# tests need no Stan; three gated tests (log-likelihood vs pbeta/dbeta, refusal
# of a mismatched cached fit, z/CTmax recovery) use a cached censored fixture.

raw_cens <- function() {
  data.frame(
    temp_C = rep(c(30, 32, 34), each = 5),
    hours  = rep(c(0.5, 1, 2, 4, 8), times = 3),
    fvfm   = c(0.95, 0.90, 0.70, 0.20, 0.00,
               0.92, 0.60, 0.30, 0.00, 0.00,
               0.80, 0.40, 0.10, 0.00, 1.00),
    lod    = c(rep(0.1, 10), rep(NA, 5)),
    # per-row limits that differ by row
    lod_var = c(0.05, 0.05, 0.05, 0.25, 0.05, 0.05, 0.05, 0.05, 0.10, 0.05,
                0.05, 0.05, 0.05, 0.05, 0.05),
    hi     = c(0.94, 0.99, 0.99, 0.99, 0.99, 0.91, 0.99, 0.99, 0.99, 0.99,
               0.99, 0.99, 0.99, 0.99, 0.99)
  )
}
std_cens <- function(limit, ...) {
  standardize_data(raw_cens(), temp = "temp_C", duration = "hours",
                   proportion = "fvfm", detection_limit = limit, ...)
}

test_that("a lower limit flags readings below it as left-censored at the limit", {
  d <- std_cens(0.05)
  zero <- d$fvfm == 0
  expect_equal(d$censored, ifelse(zero, "left", "none"))
  expect_equal(d$survival[zero], rep(0.05, sum(zero)))
  expect_equal(d$survival_recorded, raw_cens()$fvfm)
  # Uncensored rows are exactly what the uncensored path gives (clamped as before).
  plain <- standardize_data(raw_cens(), temp = "temp_C", duration = "hours",
                            proportion = "fvfm")
  expect_equal(d$survival[!zero], plain$survival[!zero])
  expect_false(any(c("censored", "survival_recorded") %in% names(plain)))
  meta <- attr(d, "tdt_meta")
  expect_true(meta$censored)
  expect_equal(meta$detection_limit, 0.05)
  expect_equal(meta$n_censored, c(left = sum(zero), right = 0))
  expect_false(attr(plain, "tdt_meta")$censored)
  expect_null(attr(plain, "tdt_meta")$n_censored)
})

test_that("a reading at the limit is observed, also after floating-point error", {
  raw <- raw_cens()
  raw$fvfm[3] <- 0.081 / 0.81           # 0.09999999999999999 in floating point
  d <- standardize_data(raw, temp = "temp_C", duration = "hours",
                        proportion = "fvfm", detection_limit = 0.1)
  expect_equal(d$censored[3], "none")
  expect_equal(d$survival[3], raw$fvfm[3])
})

test_that("positive readings below a lower limit are coarsened, with a warning", {
  expect_warning(d <- std_cens(0.25), "2 positive reading")
  expect_equal(d$censored[d$fvfm == 0.20], "left")
  expect_equal(d$censored[d$fvfm == 0.10], "left")
  expect_equal(d$survival[d$fvfm == 0.20], 0.25)
})

test_that("an upper limit right-censors readings at or above it", {
  d <- std_cens(c(lower = 0.05, upper = 0.99))
  one <- d$fvfm == 1
  expect_equal(d$censored[one], "right")
  expect_equal(d$survival[one], 0.99)
  expect_equal(attr(d, "tdt_meta")$n_censored, c(left = 4, right = 1))
  # Positional c(lower, upper) means the same.
  expect_equal(std_cens(c(0.05, 0.99))$censored, d$censored)
  # A reading exactly at the upper limit is censored ("at or above").
  raw <- raw_cens()
  raw$fvfm[1] <- 0.99
  expect_no_warning(d_at <- standardize_data(raw, temp = "temp_C", duration = "hours",
                                             proportion = "fvfm",
                                             detection_limit = c(NA, 0.99)))
  expect_equal(d_at$censored[1], "right")
  # Readings below 1 at or above the upper limit are coarsened, with a warning.
  expect_warning(d2 <- std_cens(c(NA, 0.9)), "2 reading\\(s\\) lie above the upper")
  expect_equal(d2$censored[d2$fvfm == 0.92], "right")
  expect_equal(d2$censored[d2$fvfm == 0.95], "right")
})

test_that("an upper-only limit leaves zeros clamped, not censored", {
  d <- std_cens(c(upper = 0.99))
  expect_equal(d$censored[d$fvfm == 0], rep("none", sum(d$fvfm == 0)))
  expect_equal(d$survival[d$fvfm == 0], rep(0.001, sum(d$fvfm == 0)))
  expect_equal(d$censored[d$fvfm == 1], "right")
})

test_that("per-row limits come from a column; NA means no limit for that row", {
  d <- std_cens("lod")
  expect_equal(d$censored[1:10], ifelse(raw_cens()$fvfm[1:10] == 0, "left", "none"))
  # Rows 11-15 have no limit: their zero is clamped, not censored.
  expect_equal(d$censored[11:15], rep("none", 5))
  expect_equal(d$survival[14], 0.001)
  # A list mixes a column (lower) with a number (upper).
  d2 <- std_cens(list(lower = "lod", upper = 0.99))
  expect_equal(d2$censored[15], "right")
  expect_equal(d2$censored[1:10], d$censored[1:10])
})

test_that("detection_limit input errors are clear", {
  raw <- raw_cens()
  expect_error(standardize_data(transform(raw, n = 30L, alive = 10L),
                                temp = "temp_C", duration = "hours",
                                n_total = "n", n_surv = "alive",
                                detection_limit = 0.1),
               "proportion` response only")
  expect_error(std_cens(0), "strictly between 0 and 1")
  expect_error(std_cens(1.2), "strictly between 0 and 1")
  expect_error(std_cens(c(0.5, 0.2)), "below the upper one")
  expect_error(std_cens(c(low = 0.1)), "lower` and/or `upper")
  expect_error(std_cens(c(NA_real_, NA_real_)), "sets no limit")
  expect_error(std_cens("not_a_column"), "not_a_column")
  expect_error(std_cens(c(0.1, 0.5, 0.9)), "length 1")
  expect_error(std_cens(list(lower = c(0.1, 0.2))), "one number")
})

test_that("an existing `censored` column is overwritten with a warning", {
  raw <- raw_cens()
  raw$censored <- "something else"
  expect_warning(standardize_data(raw, temp = "temp_C", duration = "hours",
                                  proportion = "fvfm", detection_limit = 0.05),
                 "overwrites existing column\\(s\\).*censored")
  raw2 <- raw_cens()
  raw2$survival_recorded <- 1
  expect_warning(standardize_data(raw2, temp = "temp_C", duration = "hours",
                                  proportion = "fvfm", detection_limit = 0.05),
                 "overwrites existing column\\(s\\).*survival_recorded")
})

test_that("the formula gets cens() on the response side for every curve and parameterisation", {
  beta <- brms::Beta(link = "identity")
  for (m in names(tdt_models())) {
    for (direct in c(FALSE, TRUE)) {
      args <- list(model = m, family = beta)
      if (direct) args <- c(args, list(ctmax = ~ 1, z = ~ 1))
      plain <- do.call(make_tdt_formula, args)
      cens  <- do.call(make_tdt_formula, c(args, list(censored = TRUE)))
      lhs   <- function(f) paste(deparse(f$formula[[2]]), collapse = "")
      rhs   <- function(f) paste(deparse(f$formula[[3]]), collapse = "")
      expect_equal(lhs(plain), "survival")
      expect_equal(lhs(cens), "survival | cens(censored)")
      expect_identical(rhs(cens), rhs(plain))
      expect_identical(lapply(cens$pforms, deparse), lapply(plain$pforms, deparse))
      expect_identical(tdt_model_from_formula(cens), m)
    }
  }
  expect_error(make_tdt_formula(censored = TRUE), "needs the Beta family")
})

test_that("fit_4pl() picks up censoring from the data and the Stan code uses the CDF", {
  d  <- std_cens(c(0.05, 0.99))
  wf <- fit_4pl(d, temp_effects = "mid", fit = FALSE)
  expect_true(wf$meta$censored)
  expect_equal(paste(deparse(wf$formula$formula[[2]]), collapse = ""),
               "survival | cens(censored)")
  code <- suppressWarnings(brms::make_stancode(wf$formula, data = wf$data,
                                               prior = wf$prior))
  expect_match(code, "beta_lcdf(", fixed = TRUE)     # below the lower limit
  expect_match(code, "beta_lccdf(", fixed = TRUE)    # above the upper limit
  # The uncensored data give the plain likelihood.
  plain <- fit_4pl(standardize_data(raw_cens(), temp = "temp_C", duration = "hours",
                                    proportion = "fvfm"),
                   temp_effects = "mid", fit = FALSE)
  expect_false(plain$meta$censored)
  expect_false(grepl("cens", paste(deparse(plain$formula$formula), collapse = "")))
})

test_that("print() reports the censored counts", {
  wf <- fit_4pl(std_cens(c(0.05, 0.99)), temp_effects = "mid", fit = FALSE)
  expect_output(print(wf), "Censored: 4 below and 1 above the detection limit")
  plain <- fit_4pl(standardize_data(raw_cens(), temp = "temp_C", duration = "hours",
                                    proportion = "fvfm"), temp_effects = "mid", fit = FALSE)
  expect_false(any(grepl("Censored", utils::capture.output(print(plain)))))
})

test_that("censored data refuse a count family and a missing flag column", {
  d <- std_cens(0.05)
  expect_error(fit_4pl(d, family = brms::beta_binomial(link = "identity"),
                       fit = FALSE), "need the Beta family")
  d2 <- d
  d2$censored <- NULL
  attr(d2, "tdt_meta") <- attr(d, "tdt_meta")
  expect_error(fit_4pl(d2, fit = FALSE), "`censored` column is missing")
})

test_that("the tls() grid marks a censored fit's rows uncensored, whatever row 1 holds", {
  d <- std_cens(0.05)
  d <- d[order(d$censored == "none"), ]         # a censored row first
  expect_equal(d$censored[1], "left")
  g <- tls_build_grid(as.data.frame(d), censored = TRUE)
  expect_true(all(g$censored == "none"))
  # An uncensored fit's own `censored` column, or a moderator of that name,
  # is left alone.
  expect_true(all(tls_build_grid(as.data.frame(d))$censored == "left"))
  m <- data.frame(temp_c = c(-1, 1), censored = c("ctrl", "heat"))
  gm <- tls_build_grid(m, by = "censored", temp_grid = 0, censored = TRUE)
  expect_setequal(gm$censored, c("ctrl", "heat"))
})

test_that("observed overlays and summaries show readings as recorded", {
  d <- std_cens(0.05)
  pred <- list(summary = data.frame(temp = 30, duration = c(0.5, 8),
                                    survival_lower = 0.1, survival_median = 0.5,
                                    survival_upper = 0.9))
  p  <- plot_survival_curves(pred, observed = d, clip_to_observed = FALSE)
  pts <- ggplot2::layer_data(p, 3)
  expect_equal(sort(pts$y), sort(d$survival_recorded))
  s <- summarise_observed_survival(d)
  cell <- s[s$temp == 32 & s$duration == 8, ]
  expect_equal(cell$survival_mean, 0)            # the censored 0, not the 0.05 limit
})

test_that("a censored fit's log-likelihood is the Beta CDF at the limit (gated)", {
  skip_unless_brms()
  wf  <- load_fixture_workflow_beta_censored()
  fit <- get_brmsfit(wf)
  ids <- 1:50
  ll  <- brms::log_lik(fit, draw_ids = ids)
  mu  <- brms::posterior_linpred(fit, transform = TRUE, draw_ids = ids)
  phi <- posterior::as_draws_df(fit)$phi[ids]
  y   <- wf$data$survival
  cen <- wf$data$censored == "left"
  expect_gt(sum(cen), 10)
  a <- mu * phi
  b <- (1 - mu) * phi
  hand <- matrix(stats::dbeta(rep(y, each = length(ids)), a, b, log = TRUE),
                 nrow = length(ids))
  hand[, cen] <- stats::pbeta(rep(y[cen], each = length(ids)), a[, cen], b[, cen],
                              log.p = TRUE)
  expect_equal(unname(ll), unname(hand), tolerance = 1e-8)
})

test_that("a cached censored fit is refused when uncensored readings are requested (gated)", {
  skip_unless_brms()
  load_fixture_workflow_beta_censored()           # makes sure the cached file exists
  raw <- simulate_tdt_beta(durations_hours = phase3_hours, n_rep = 8, floor = 0.1,
                           seed = 20261126)
  plain <- standardize_data(raw, temp = "T", duration = "t_hours",
                            proportion = "y_prop", duration_unit = "hours")
  expect_error(
    fit_4pl(plain, chains = 2, iter = 1500, warmup = 750, cores = 2, seed = 1,
            file = here::here("tests", "testthat", "fixtures", "sim_fit_beta_censored"),
            file_refit = "never", refresh = 0),
    "has censored readings")
})

test_that("a censored fit recovers z and CTmax when about half the readings sit below the floor (gated)", {
  skip_unless_brms()
  wf <- load_fixture_workflow_beta_censored()
  tr <- attr(wf, "truth")
  expect_gt(sum(wf$data$censored == "left"), 50)
  s <- suppressWarnings(tls(wf, params = c("z", "ctmax"), t_ref = 60, lethal = FALSE,
                            probs = c(0.001, 0.5, 0.999)))$summary
  g <- function(q) s[s$quantity == q, ]
  # Tolerances from ten seeds of this design (censoring_sim in the notes:
  # largest errors z 0.89, CTmax 0.34 with 4 x 1000 draws).
  expect_lt(abs(g("z")$median - tr$z), 1.5)
  expect_lt(abs(g("CTmax")$median - tr$CTmax_rel_1hr), 0.6)
  expect_true(g("z")$lower <= tr$z && tr$z <= g("z")$upper)
  expect_true(g("CTmax")$lower <= tr$CTmax_rel_1hr && tr$CTmax_rel_1hr <= g("CTmax")$upper)
  # R^2 warns (once, in the package's words) that it is only a rough guide.
  w <- testthat::capture_warnings(bayes_R2_tls(wf))
  expect_length(w, 1)
  expect_match(w, "rough guide")
})


# --- Review fixes (2026-11-26) ------------------------------------------------

test_that("an NA (of any type) on one side means no limit there; empty names error", {
  ref_lo <- std_cens("lod")
  expect_equal(std_cens(c(lower = "lod", upper = NA))$censored, ref_lo$censored)
  expect_equal(nrow(std_cens(c(lower = "lod", upper = NA))), 15)
  expect_warning(ref_hi <- std_cens(c(upper = "hi")), "above the upper")
  expect_warning(d <- std_cens(c(NA, "hi")), "above the upper")
  expect_equal(d$censored, ref_hi$censored)
  expect_equal(std_cens(list(lower = 0.1, upper = NA_character_))$censored,
               std_cens(0.1)$censored)
  expect_error(std_cens(c(lower = "lod", upper = "")), "not in the data")
})

test_that("per-row limits that differ by row are applied row by row", {
  expect_warning(d <- std_cens("lod_var"), "positive reading")
  left <- d$censored == "left"
  expect_equal(d$survival[left], raw_cens()$lod_var[left])
  expect_equal(d$censored[4], "left")             # 0.20 below its own 0.25 limit
  expect_equal(d$censored[13], "none")            # 0.10 above its 0.05 limit
  # A per-row upper column, alone or mixed with a number or a lower column.
  expect_warning(d2 <- std_cens(list(lower = 0.05, upper = "hi")), "2 reading\\(s\\) lie above")
  right <- d2$censored == "right"
  expect_equal(which(right), c(1, 6, 15))
  expect_equal(d2$survival[right], raw_cens()$hi[right])
  expect_warning(d3 <- std_cens(c(lower = "lod", upper = "hi")), "lie above")
  expect_equal(attr(d3, "tdt_meta")$n_censored, c(left = 3, right = 3))
})

test_that("limit columns must be numeric", {
  raw <- raw_cens()
  raw$lod_chr <- ifelse(is.na(raw$lod), NA, "0,1")
  expect_error(standardize_data(raw, temp = "temp_C", duration = "hours",
                                proportion = "fvfm", detection_limit = "lod_chr"),
               "must be numeric")
  raw$lod_fac <- factor(raw$lod)
  expect_error(standardize_data(raw, temp = "temp_C", duration = "hours",
                                proportion = "fvfm", detection_limit = "lod_fac"),
               "must be numeric")
})

test_that("proportion_eps must be one number in (0, 0.5)", {
  for (eps in list(NULL, NA, 0, 0.6, c(0.01, 0.02))) {
    expect_error(std_cens(0.05, proportion_eps = eps), "proportion_eps")
  }
  expect_error(standardize_data(raw_cens(), temp = "temp_C", duration = "hours",
                                proportion = "fvfm", proportion_eps = NULL),
               "proportion_eps")
})

test_that("NA readings and rows dropped in cleaning keep everything aligned", {
  raw <- raw_cens()
  raw$fvfm[2] <- NA
  raw$hours[10] <- 0                                  # a censored zero, dropped
  raw$temp_C[14] <- NA                                # another censored zero, dropped
  raw$fvfm[11] <- 0.02; raw$hours[11] <- 0            # positive below limit, dropped
  expect_no_warning(d <- standardize_data(raw, temp = "temp_C", duration = "hours",
                                          proportion = "fvfm", detection_limit = 0.05))
  expect_equal(nrow(d), 11)
  expect_false(anyNA(d$censored))
  expect_equal(attr(d, "tdt_meta")$n_censored, c(left = 2, right = 0))
})

test_that("limits are stored as given and only unlimited sides are clamped", {
  raw <- raw_cens()
  raw$fvfm[3] <- 0.07
  d <- standardize_data(raw, temp = "temp_C", duration = "hours", proportion = "fvfm",
                        detection_limit = 0.05, proportion_eps = 0.1)
  expect_equal(d$survival[3], 0.07)                   # observed, above the limit: kept
  expect_equal(d$survival[15], 0.9)                   # the 1 has no upper limit: clamped
  d2 <- std_cens(c(0.0005, 0.9995))
  expect_equal(unique(d2$survival[d2$censored == "left"]), 0.0005)
  expect_equal(unique(d2$survival[d2$censored == "right"]), 0.9995)
  d3 <- std_cens(1e-8)                                # a tiny limit still censors zeros
  expect_equal(unique(d3$survival[d3$fvfm == 0]), 1e-8)
  expect_true(all(d3$censored[d3$fvfm == 0] == "left"))
  raw$fvfm[1] <- 0.8019 / 0.81                        # 0.98999999999999988: at 0.99
  expect_no_warning(d4 <- standardize_data(raw, temp = "temp_C", duration = "hours",
                                           proportion = "fvfm",
                                           detection_limit = c(0.05, 0.99)))
  expect_equal(d4$censored[1], "right")
})

test_that("new arguments come last, so positional calls keep their meaning", {
  expect_equal(tail(names(formals(standardize_data)), 1), "detection_limit")
  expect_equal(tail(names(formals(make_4pl_formula)), 1), "censored")
  expect_equal(tail(names(formals(make_tdt_formula)), 1), "censored")
  f <- make_4pl_formula(NULL, NULL, brms::beta_binomial(link = "identity"),
                        "survival", "mid")
  expect_equal(paste(deparse(f$pforms$lowraw), collapse = ""), "lowraw ~ 1")
  expect_error(make_tdt_formula(family = brms::Beta(link = "identity"),
                                censored = "yes"), "TRUE or FALSE")
  expect_error(make_tdt_formula(family = brms::Beta(link = "identity"),
                                censored = NA), "TRUE or FALSE")
})

test_that("censoring flags without their metadata are refused", {
  d <- std_cens(0.05)
  m <- merge(d, data.frame(temp = c(30, 32, 34), site = "A"))
  expect_null(attr(m, "tdt_meta"))
  expect_error(fit_4pl(m, family = brms::Beta(link = "identity"),
                       temp_effects = "mid", fit = FALSE),
               "attribute was lost")
})

test_that("a cached fit whose censoring disagrees with the request is refused", {
  fc <- brms::bf(survival | cens(censored) ~ a, a ~ 1, nl = TRUE)
  fu <- brms::bf(survival ~ a, a ~ 1, nl = TRUE)
  expect_true(tdt_check_cached_censoring(fc, TRUE))
  expect_true(tdt_check_cached_censoring(fu, FALSE))
  expect_error(tdt_check_cached_censoring(fu, TRUE, "x"), "has uncensored readings")
  expect_error(tdt_check_cached_censoring(fc, FALSE, "x"), "has censored readings")
})

test_that("every fitter and the logit link pick up censoring", {
  d <- std_cens(c(0.05, 0.99))
  lhs <- function(wf) paste(deparse(wf$formula$formula[[2]]), collapse = "")
  for (wf in list(fit_lognormal(d, temp_effects = "mid", fit = FALSE),
                  fit_weibull(d, type = 1, temp_effects = "mid", fit = FALSE),
                  fit_weibull(d, type = 2, temp_effects = "mid", fit = FALSE),
                  fit_4pl(d, family = brms::Beta(link = "logit"),
                          temp_effects = "mid", fit = FALSE))) {
    expect_true(wf$meta$censored)
    expect_equal(lhs(wf), "survival | cens(censored)")
  }
})

test_that("brms receives the censoring codes and the limits", {
  d  <- std_cens(c(0.05, 0.99))
  wf <- fit_4pl(d, temp_effects = "mid", fit = FALSE)
  sd <- suppressWarnings(brms::make_standata(wf$formula, data = wf$data, prior = wf$prior))
  expect_equal(as.vector(sd$cens), unname(c(none = 0, left = -1, right = 1)[d$censored]))
  expect_equal(as.vector(sd$Y)[d$censored != "none"],
               d$survival[d$censored != "none"])
  expect_equal(sort(unique(as.vector(sd$Y)[d$censored != "none"])), c(0.05, 0.99))
})

test_that("print() counts the censored rows actually in the workflow", {
  d <- std_cens(c(0.05, 0.99))
  sub <- d[d$temp != 34, ]                             # keeps tdt_meta
  expect_equal(attr(sub, "tdt_meta")$n_censored, c(left = 4, right = 1))
  wf <- fit_4pl(sub, temp_effects = "mid", fit = FALSE)
  expect_output(print(wf), "Censored: 3 below and 0 above")
  expect_equal(wf$meta$n_censored, c(left = 3, right = 0))
})

test_that("bayes_R2_tls() warns once, and only when readings are censored", {
  fake_r2 <- function(object, ...) {
    warning("Results may not be meaningful for censored models.")
    matrix(c(0.8, 0.01, 0.78, 0.82), 1,
           dimnames = list(NULL, c("Estimate", "Est.Error", "Q2.5", "Q97.5")))
  }
  testthat::local_mocked_bindings(bayes_R2 = fake_r2, .package = "brms")
  wf <- fit_4pl(std_cens(0.05), temp_effects = "mid", fit = FALSE)
  wf$fit <- structure(list(), class = "fake_fit")
  w <- testthat::capture_warnings(r <- bayes_R2_tls(wf))
  expect_length(w, 1)
  expect_match(w, "rough guide")
  expect_equal(r$estimate, 0.8)
  # Censoring switched on but no row censored: no package warning (brms's
  # generic note is still muffled).
  raw <- raw_cens()
  raw$fvfm[raw$fvfm == 0] <- 0.08
  wf0 <- fit_4pl(standardize_data(raw, temp = "temp_C", duration = "hours",
                                  proportion = "fvfm", detection_limit = 0.05),
                 temp_effects = "mid", fit = FALSE)
  wf0$fit <- wf$fit
  expect_length(testthat::capture_warnings(bayes_R2_tls(wf0)), 0)
  # Uncensored fit: brms's own warning (if any) passes through untouched.
  wfu <- fit_4pl(standardize_data(raw_cens(), temp = "temp_C", duration = "hours",
                                  proportion = "fvfm"), temp_effects = "mid", fit = FALSE)
  wfu$fit <- wf$fit
  expect_match(testthat::capture_warnings(bayes_R2_tls(wfu)), "meaningful")
})

test_that("the prediction grid and the landscape overlay handle censored data", {
  d  <- std_cens(0.05)
  wf <- fit_4pl(d, temp_effects = "mid", fit = FALSE)
  g  <- new_tdt_grid(wf, temps = c(30, 34), durations = c(1, 2))
  expect_true("censored" %in% names(g))
  expect_true(all(g$censored == "none"))
  grid <- expand.grid(temp = seq(30, 38, length.out = 12),
                      duration = 10^seq(log10(0.5), log10(30), length.out = 12))
  grid$survival_median <- plogis(8 - 0.4 * grid$duration - 0.5 * (grid$temp - 30))
  p <- plot_tdt_landscape(list(summary = tibble::as_tibble(grid)), observed = d,
                          clip_to_observed = FALSE)
  pts <- ggplot2::ggplot_build(p)$data
  pts <- pts[[length(pts)]]
  expect_equal(nrow(pts), nrow(d))
  # The fill is the reading as recorded (0), not the 0.05 limit.
  ref <- d
  ref$survival <- ref$survival_recorded
  ref$censored <- NULL
  p_ref <- plot_tdt_landscape(list(summary = tibble::as_tibble(grid)), observed = ref,
                              clip_to_observed = FALSE)
  pts_ref <- ggplot2::ggplot_build(p_ref)$data
  expect_equal(pts$fill, pts_ref[[length(pts_ref)]]$fill)
  # An unrelated raw `survival_recorded` column on uncensored data is ignored.
  plain <- standardize_data(raw_cens(), temp = "temp_C", duration = "hours",
                            proportion = "fvfm")
  plain$survival_recorded <- 1
  s <- summarise_observed_survival(plain)
  expect_equal(s$survival_mean[s$temp == 32 & s$duration == 8], 0.001)
})

# --- Second review round (2026-11-26) ----------------------------------------

test_that("a user's own `censored` column does not block an uncensored fit", {
  raw <- data.frame(temp_C = rep(c(30, 32, 34), each = 4), hours = rep(c(1, 2, 4, 8), 3),
                    n = 20L, alive = c(20, 18, 12, 5, 19, 14, 6, 1, 17, 8, 2, 0),
                    censored = rep(c("none", "none", "none", "right"), 3))
  std <- standardize_data(raw, "temp_C", "hours", n_total = "n", n_surv = "alive")
  wf  <- fit_4pl(std, temp_effects = "mid", fit = FALSE)
  expect_false(wf$meta$censored)
  expect_equal(paste(deparse(wf$formula$formula[[2]]), collapse = ""),
               "n_surv | trials(n_total)")
  # A grouped uncensored fit whose moderator is called `censored` keeps it.
  raw$grp <- rep(c("ctrl", "heat"), 6)
  names(raw)[names(raw) == "censored"] <- "code"
  names(raw)[names(raw) == "grp"] <- "censored"
  std2 <- standardize_data(raw, "temp_C", "hours", n_total = "n", n_surv = "alive")
  wf2  <- fit_4pl(std2, by = "censored", fit = FALSE)
  g <- new_tdt_grid(wf2, temps = 30, durations = 1, by = "censored")
  expect_setequal(as.character(g$censored), c("ctrl", "heat"))
})

test_that("data re-standardised without a limit keep a plain likelihood and recorded summaries", {
  d <- std_cens(0.05)
  plain <- suppressWarnings(standardize_data(d, temp = "temp", duration = "duration",
                                             proportion = "survival_recorded"))
  expect_false(attr(plain, "tdt_meta")$censored)
  wf <- fit_4pl(plain, temp_effects = "mid", fit = FALSE)
  expect_equal(paste(deparse(wf$formula$formula[[2]]), collapse = ""), "survival")
  s <- summarise_observed_survival(plain)
  expect_equal(s$survival_mean[s$temp == 32 & s$duration == 8], 0.001)  # clamped, not swapped
})

test_that("binding censored and uncensored standardised data is refused, not silently thinned", {
  a <- std_cens(0.05, temp_mean = 32)
  b <- standardize_data(raw_cens(), temp = "temp_C", duration = "hours",
                        proportion = "fvfm", temp_mean = 32)
  ab <- dplyr::bind_rows(a, b)
  expect_true(anyNA(ab$censored))
  expect_error(fit_4pl(ab, temp_effects = "mid", fit = FALSE), "other than")
  # The observed overlay keeps `survival` where nothing was recorded separately.
  o <- tdt_observed_as_recorded(ab)
  expect_false(anyNA(o$survival))
})

test_that("floating-point noise at 0 and 1 raises no coarsening warning", {
  raw <- raw_cens()
  raw$fvfm[15] <- 0.3 / (0.1 * 3)                     # 0.99999999999999978
  expect_no_warning(d <- standardize_data(raw, temp = "temp_C", duration = "hours",
                                          proportion = "fvfm",
                                          detection_limit = c(0.05, 0.99)))
  expect_equal(d$censored[15], "right")
  raw <- raw_cens()
  raw$fvfm[5] <- 1 - 0.3 / (0.1 * 3)                  # 2.2e-16
  expect_no_warning(d2 <- standardize_data(raw, temp = "temp_C", duration = "hours",
                                           proportion = "fvfm", detection_limit = 0.05))
  expect_equal(d2$censored[5], "left")
})

test_that("zero rows with a valid limit reach the usual empty-data warning", {
  raw <- raw_cens()[0, ]
  expect_warning(d <- standardize_data(raw, temp = "temp_C", duration = "hours",
                                       proportion = "fvfm", detection_limit = 0.1),
                 "Only 0 unique")
  expect_equal(nrow(d), 0)
})
