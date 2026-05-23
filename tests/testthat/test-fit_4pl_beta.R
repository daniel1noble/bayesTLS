# Gated brms integration test for the Beta (continuous-proportion) path:
# standardize_data(proportion =) -> fit_4pl(family = Beta(link = "identity")).
# Asserts z recovery within tolerance of truth and that the downstream helpers
# that used to assume an n_total column now run on a Beta fit. Set
# RUN_BRMS_TESTS=true to enable; otherwise a no-op.

test_that("fit_4pl recovers a known Beta truth and records the family", {
  skip_unless_brms()

  wf <- load_fixture_workflow_beta()
  ts <- truth_summary()

  # The workflow should know it is a proportion / Beta fit with no n_total.
  expect_equal(wf$meta$response_type, "proportion")
  expect_equal(wf$meta$family, "beta")
  expect_equal(wf$meta$link,   "identity")
  expect_false("n_total" %in% names(wf$data))

  out <- extract_tdt(wf, t_ref = 60, ndraws = 500, lethal = FALSE)

  # Posterior median within tolerance of truth and truth inside the 95% CrI.
  expect_lt(abs(out$z$summary$z_median        - ts$z),         1.5)
  expect_lt(abs(out$CTmax$summary$temp_median - ts$CTmax_1hr), 1.5)
  expect_lte(out$z$summary$z_lower, ts$z)
  expect_gte(out$z$summary$z_upper, ts$z)

  # T_crit only makes sense for lethal endpoints; off by default here.
  expect_null(out$T_crit)
  expect_false(out$meta$lethal)
})

test_that("downstream helpers run on a Beta fit (no n_total)", {
  skip_unless_brms()

  wf <- load_fixture_workflow_beta()

  # predict_survival_curves: finite proportions in (0, 1).
  psc <- predict_survival_curves(wf, temps = c(30, 33, 36), ndraws = 200)
  expect_true(all(is.finite(psc$summary$survival_median)))
  expect_true(all(psc$summary$survival_median > 0 &
                  psc$summary$survival_median < 1))
  # The grid must NOT carry an n_total column for a Beta fit.
  expect_false("n_total" %in% names(psc$grid))

  # derive_tdt_curve: time-to-threshold decreases with temperature.
  tc <- derive_tdt_curve(wf, temp_grid = c(30, 33, 36), ndraws = 200,
                         time_multiplier = 60, output_time_unit = "min")
  expect_true(all(diff(tc$summary$duration_median) < 0))

  # derive_tdt_landscape: finite survival surface, no NAs.
  ls <- derive_tdt_landscape(wf,
                             temp_grid     = seq(30, 36, length.out = 15),
                             duration_grid = seq(0.05, 50, length.out = 15),
                             ndraws        = 150)
  expect_false(any(is.na(ls$summary$survival_median)))

  # summarise_observed_survival: runs without an n_total column and omits the
  # n_total_sum summary.
  sos <- summarise_observed_survival(wf$data)
  expect_true(all(is.finite(sos$survival_mean)))
  expect_false("n_total_sum" %in% names(sos))
})
