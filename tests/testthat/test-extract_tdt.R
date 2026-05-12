# Gated brms integration test: fits a small cached model and asserts that
# extract_tdt() recovers z, CTmax_1hr, and T_crit within tolerance of truth.
# Set RUN_BRMS_TESTS=true to enable. Without that, the file becomes a no-op.

test_that("extract_tdt recovers z, CTmax_1hr, and T_crit from simulated truth", {
  skip_unless_brms()

  wf      <- load_fixture_workflow()
  ts      <- truth_summary()
  out     <- extract_tdt(wf, t_ref = 60, TC_thresh = 0.05, ndraws = 500)

  # Posterior medians within 1 °C of truth — tight enough to detect bugs,
  # loose enough that random MCMC noise on a small fit won't trigger false
  # failures.
  expect_lt(abs(out$z$summary$z_median             - ts$z),         1.5)
  expect_lt(abs(out$CTmax$summary$temp_median      - ts$CTmax_1hr), 1.0)

  # Truth should sit inside the 95% credible intervals.
  expect_lte(out$z$summary$z_lower,         ts$z)
  expect_gte(out$z$summary$z_upper,         ts$z)
  expect_lte(out$CTmax$summary$temp_lower,  ts$CTmax_1hr)
  expect_gte(out$CTmax$summary$temp_upper,  ts$CTmax_1hr)

  # T_crit should always sit below CTmax_1hr (LD5 temperature is lower than
  # LD50 temperature at the same exposure duration).
  expect_lt(out$T_crit$summary$temp_median, out$CTmax$summary$temp_median)
})

test_that("predict_heat_injury recovers analytical planted dose within CrI", {
  skip_unless_brms()

  wf <- load_fixture_workflow()
  ts <- truth_summary()

  scens <- make_temperature_scenarios(baseline = 20, spike_temp = 28,
                                      n_hours = 96,
                                      spike_times_single = 24,
                                      spike_times_multi  = c(24, 48, 72))
  T_c_val <- 24

  planted_single <- planted_dose_from_trace(
    scens$single_spike, z = ts$z, CTmax_1hr = ts$CTmax_1hr, T_c = T_c_val
  )
  planted_multi  <- planted_dose_from_trace(
    scens$multi_spike, z = ts$z, CTmax_1hr = ts$CTmax_1hr, T_c = T_c_val
  )

  hi_single <- predict_heat_injury(scens$single_spike, wf,
                                   T_c = T_c_val, ndraws = 300)
  hi_multi  <- predict_heat_injury(scens$multi_spike,  wf,
                                   T_c = T_c_val, ndraws = 300)

  truth_single <- tail(planted_single$hi_cumulative, 1)
  truth_multi  <- tail(planted_multi$hi_cumulative,  1)
  final_single <- tail(hi_single$summary, 1)
  final_multi  <- tail(hi_multi$summary,  1)

  # Posterior median within ~30% relative (small fit + ~10 spikes); truth
  # inside 95% CrI for both scenarios.
  expect_lt(abs(final_single$hi_median - truth_single) / truth_single, 0.3)
  expect_lt(abs(final_multi$hi_median  - truth_multi)  / truth_multi,  0.3)

  expect_lte(final_single$hi_lower, truth_single)
  expect_gte(final_single$hi_upper, truth_single)
  expect_lte(final_multi$hi_lower,  truth_multi)
  expect_gte(final_multi$hi_upper,  truth_multi)

  # Flat trace at sub-T_c baseline → zero HI by construction.
  hi_flat <- predict_heat_injury(scens$flat, wf,
                                 T_c = T_c_val, ndraws = 100)
  expect_equal(tail(hi_flat$summary$hi_median, 1), 0)
})
