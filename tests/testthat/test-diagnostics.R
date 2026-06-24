# Tests for the sampling-diagnostic helpers (R/diagnostics.R). The NULL-fit
# guards are brms-free; the substantive checks use the cached fixture fit and
# assert internal consistency (e.g. all_pass is exactly the conjunction of the
# component flags) and cross-validate tdt_parameter_table()'s z against
# extract_tdt(). Set RUN_BRMS_TESTS=true to enable the gated block.

test_that("diagnose_tdt_fit and tdt_parameter_table reject an unfitted workflow", {
  expect_error(diagnose_tdt_fit(list(fit = NULL)),    "NULL")
  expect_error(tdt_parameter_table(list(fit = NULL)), "NULL")
})

test_that("diagnose_tdt_fit returns one internally-consistent diagnostic row", {
  skip_unless_brms()
  wf <- load_fixture_workflow()
  d  <- diagnose_tdt_fit(wf)

  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 1L)
  expect_true(all(c("rhat_max", "ess_bulk_min", "ess_tail_min", "divergences",
                    "treedepth_max", "treedepth_hits", "bfmi_min", "rhat_pass",
                    "ess_pass", "divergence_pass", "treedepth_pass", "bfmi_pass",
                    "all_pass") %in% names(d)))

  # Statistic ranges that must hold for any valid fit.
  expect_true(is.finite(d$rhat_max) && d$rhat_max >= 1)   # Rhat is >= 1
  expect_gt(d$ess_bulk_min, 0)
  expect_gt(d$ess_tail_min, 0)
  expect_gte(d$divergences, 0L)
  expect_gte(d$treedepth_hits, 0L)
  # treedepth_max is the CONFIGURED ceiling read from the fit; saturations are
  # counted against it. A healthy fixture that never reaches the ceiling reports
  # zero hits -- regression guard for the old `sum(td_vals >= max(td_vals))`,
  # which always flagged >= 1 (the observed max is always attained).
  expect_gte(d$treedepth_max, 1L)
  expect_equal(d$treedepth_hits, 0L)

  # Pass flags are logical scalars and encode their stated thresholds.
  for (f in c("rhat_pass", "ess_pass", "divergence_pass",
              "treedepth_pass", "bfmi_pass", "all_pass")) {
    expect_type(d[[f]], "logical"); expect_length(d[[f]], 1L)
  }
  expect_equal(d$rhat_pass,       d$rhat_max < 1.01)
  expect_equal(d$ess_pass,        d$ess_bulk_min > 400 && d$ess_tail_min > 400)
  expect_equal(d$divergence_pass, d$divergences == 0L)
  expect_equal(d$treedepth_pass,  d$treedepth_hits == 0L)
  # all_pass is EXACTLY the conjunction of the five component flags.
  expect_equal(d$all_pass,
               d$rhat_pass && d$ess_pass && d$divergence_pass &&
                 d$treedepth_pass && d$bfmi_pass)
})

test_that("tdt_parameter_table returns natural-scale parameters that obey their constraints", {
  skip_unless_brms()
  wf <- load_fixture_workflow()
  tab <- tdt_parameter_table(wf)

  expect_equal(nrow(tab), 6L)
  expect_named(tab, c("parameter", "median", "lower", "upper"))
  # Every credible interval is ordered lower <= median <= upper.
  expect_true(all(tab$lower <= tab$median + 1e-9))
  expect_true(all(tab$median <= tab$upper + 1e-9))

  val <- function(p) tab$median[grepl(p, tab$parameter)]
  low_med <- val("^low"); up_med <- val("^up"); k_med <- val("^k ")
  # Upper asymptote sits above the lower one (truth: u = 0.95 > ell = 0.05).
  expect_gt(up_med, low_med)
  expect_gt(k_med, 0)                                   # slope k is positive
  # Asymptotes stay inside the reparameterisation bounds.
  b <- wf$meta$bounds
  expect_gte(low_med, b$low_min); expect_lte(low_med, b$low_min + b$low_w)
  expect_gte(up_med,  b$up_min);  expect_lte(up_med,  b$up_min  + b$up_w)

  # The z row must agree with extract_tdt()'s z (both are -1 / b_mid_temp_c).
  z_tab <- val("^z ")
  z_et  <- extract_tdt(wf, t_ref = 60, ndraws = 1000)$z$summary$z_median
  expect_equal(z_tab, z_et, tolerance = 0.25)
})

test_that("bayes_R2_tls rejects an unfitted workflow", {
  expect_error(bayes_R2_tls(list(fit = NULL)), "no fit")
})

test_that("bayes_R2_tls returns a tidy one-row R^2 summary", {
  skip_unless_brms()
  wf <- load_fixture_workflow()
  r  <- bayes_R2_tls(wf)

  expect_s3_class(r, "tbl_df")
  expect_equal(nrow(r), 1L)
  expect_named(r, c("estimate", "est_error", "lower", "upper"))
  # R^2 is a proportion of variance explained: inside [0, 1] and interval-ordered.
  expect_gte(r$estimate, 0); expect_lte(r$estimate, 1)
  expect_lte(r$lower, r$estimate); expect_gte(r$upper, r$estimate)
  expect_gt(r$est_error, 0)

  # Matches brms::bayes_R2() read directly (the wrapper only reshapes it).
  raw <- brms::bayes_R2(get_brmsfit(wf))
  expect_equal(r$estimate, unname(raw[1, "Estimate"]), tolerance = 1e-8)
})
