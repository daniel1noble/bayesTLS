# Accuracy tests for tls(). The brms-free block checks argument validation; the
# gated block validates numerical correctness against the exact -1/beta1
# identity and the simulation truth (z, CTmax, T_crit), plus the per-group
# machinery, threshold edge cases and ndraws clamping.

test_that("tls() validates arguments before touching a fit", {
  expect_error(tls(list(foo = 1)), "bayes_tls|brmsfit")
  expect_error(tls(list(foo = 1), params = "tcrit", lethal = FALSE), "lethal")
})

test_that("tls() relative z equals the exact -1/b_mid_temp_c identity and recovers truth", {
  skip_unless_brms()
  wf <- load_fixture_workflow()
  tr <- truth_summary()
  d  <- posterior::as_draws_df(get_brmsfit(wf))

  tl    <- tls(wf, params = "z", mode = "relative")
  z_tls <- tl$summary$median[tl$summary$quantity == "z"]

  # Exact closed form: relative z = -1 / b_mid_temp_c per draw (tls uses all draws).
  expect_equal(z_tls, stats::median(-1 / d$b_mid_temp_c), tolerance = 1e-6)
  # Recovers the known simulation truth (z = -1/m_beta1 = 5.556).
  expect_equal(z_tls, tr$z, tolerance = 0.5)
})

test_that("tls() recovers z and CTmax_1hr from simulated truth", {
  skip_unless_brms()
  wf  <- load_fixture_workflow()
  tr  <- truth_summary()
  out <- tls(wf, params = c("z", "ctmax"), t_ref = 60, ndraws = 500)$summary
  g   <- function(q, col) out[[col]][out$quantity == q]

  # Posterior medians within tolerance of truth ...
  expect_lt(abs(g("z", "median")     - tr$z),         1.5)
  expect_lt(abs(g("CTmax", "median") - tr$CTmax_1hr), 1.0)
  # ... and truth inside the 95% credible interval.
  expect_lte(g("z", "lower"),     tr$z)
  expect_gte(g("z", "upper"),     tr$z)
  expect_lte(g("CTmax", "lower"), tr$CTmax_1hr)
  expect_gte(g("CTmax", "upper"), tr$CTmax_1hr)
})

test_that("tls(lethal = TRUE) returns a sensible T_crit below CTmax", {
  skip_unless_brms()
  wf  <- load_fixture_workflow()
  tr  <- truth_summary()
  out <- suppressMessages(
    tls(wf, params = "all", t_ref = 60, ndraws = 500, lethal = TRUE))$summary
  g   <- function(q, col) out[[col]][out$quantity == q]

  # With default TC_rate_range = c(0.1, 1) the T_crit median sits near
  # CTmax_1hr - 2.5 * z (geometric mean of the log10(rate) range); the true
  # T_crit median must fall inside the posterior 95% CrI, and below CTmax.
  true_T_crit_median <- tr$CTmax_1hr - 2.5 * tr$z
  expect_lt(g("Tcrit", "median"), g("CTmax", "median"))
  expect_lte(g("Tcrit", "lower"), true_T_crit_median)
  expect_gte(g("Tcrit", "upper"), true_T_crit_median)
  # CrI width > 0.5 z (pooling over log10(r*) contributes ~1 z of spread).
  expect_gt(g("Tcrit", "upper") - g("Tcrit", "lower"), 0.5 * tr$z)
})

test_that("an absolute target_surv outside the fitted asymptotes warns clearly", {
  skip_unless_brms()
  wf <- load_fixture_workflow_grouped()   # survival spans ~[0.04, 0.97]
  # 0.1 lies below the lower asymptote for some draws -> LT undefined there:
  # a clear, reasoned warning instead of a bare "NaNs produced".
  expect_warning(
    tls(wf, by = "grp", target_surv = 0.1, ndraws = NULL,
        temp_mean = wf$meta$temp_mean),
    "outside the fitted curve's asymptotes")
  # an in-range absolute threshold does not raise that warning
  asymp_warn <- FALSE
  withCallingHandlers(
    tls(wf, by = "grp", target_surv = 0.5, ndraws = NULL,
        temp_mean = wf$meta$temp_mean),
    warning = function(w) {
      if (grepl("asymptotes", conditionMessage(w))) asymp_warn <<- TRUE
      invokeRestart("muffleWarning")
    })
  expect_false(asymp_warn)
})

test_that("derive_tdt_curve() clamps ndraws to the posterior size", {
  skip_unless_brms()
  wf <- load_fixture_workflow()                 # ~1500 posterior draws
  # Regression: relative-mode derive_tdt_curve() once passed ndraws straight to
  # brms::posterior_linpred(), which errors when ndraws exceeds the posterior
  # size. Requesting more draws than exist must clamp, not crash.
  expect_no_error(derive_tdt_curve(wf, temp_grid = c(30, 34, 38), ndraws = 5000))
})

test_that("tls() params/lethal switches and summary shape", {
  skip_unless_brms()
  wf <- load_fixture_workflow()

  z_only <- tls(wf, params = "z", mode = "relative")
  expect_s3_class(z_only, "tls")
  expect_equal(unique(z_only$summary$quantity), "z")
  expect_named(z_only$summary, c("quantity", "median", "lower", "upper"))

  all3 <- tls(wf, params = "all", lethal = TRUE, mode = "relative")
  expect_setequal(unique(all3$summary$quantity), c("z", "CTmax", "Tcrit"))
  # draws slot carries one row per draw per quantity
  expect_true(all(c("quantity", ".draw", "value") %in% names(all3$draws)))
})

test_that("tls() warns and NAs the summary when a group's z/CTmax draws are all non-finite", {
  skip_unless_brms()
  wf  <- load_fixture_workflow()
  fit <- get_brmsfit(wf)
  # A single-temperature grid makes the LS slope 0/0 (NaN) -> z = -1/NaN non-finite
  # for every draw, so the group's summary must be NA and a naming warning fired.
  nd  <- data.frame(grp = factor("solo"),
                    temp_c = mean(fit$data$temp_c))

  expect_warning(
    out <- tls(wf, by = "grp", newdata = nd, params = "z", mode = "relative"),
    "non-finite.*solo|solo.*non-finite"
  )
  expect_true(all(is.na(out$summary$median[out$summary$quantity == "z"])))
  # Raw draws are still returned (non-finite, not silently dropped).
  expect_true(all(!is.finite(out$draws$value[out$draws$quantity == "z"])))
})

test_that("tls() per-group machinery: a moderator the model ignores gives identical z", {
  skip_unless_brms()
  wf  <- load_fixture_workflow()
  fit <- get_brmsfit(wf)
  tc  <- seq(min(fit$data$temp_c), max(fit$data$temp_c), length.out = 5)
  # `grp` is NOT in the model, so both groups must yield the same dose-response.
  nd  <- expand.grid(grp = factor(c("a", "b")), temp_c = tc)

  byg <- tls(wf, by = "grp", newdata = nd, params = "z", mode = "relative")$summary
  za  <- byg$median[byg$grp == "a" & byg$quantity == "z"]
  zb  <- byg$median[byg$grp == "b" & byg$quantity == "z"]
  pooled <- tls(wf, params = "z", mode = "relative")$summary$median

  expect_equal(za, zb, tolerance = 1e-8)       # model ignores grp -> identical
  expect_equal(za, pooled, tolerance = 0.02)   # and equals the pooled estimate
})
