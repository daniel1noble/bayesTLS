# Gated brms integration test: fits a small cached model and asserts that
# extract_tdt() recovers z and CTmax_1hr within tolerance of truth, and that
# `lethal = TRUE` returns a sensible T_crit. Set RUN_BRMS_TESTS=true to
# enable. Without that, the file becomes a no-op.

test_that("extract_tdt recovers z and CTmax_1hr from simulated truth", {
  skip_unless_brms()

  wf      <- load_fixture_workflow()
  ts      <- truth_summary()
  out     <- extract_tdt(wf, t_ref = 60, ndraws = 500)

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

  # T_crit is opt-in via `lethal = TRUE` and must be NULL otherwise.
  expect_null(out$T_crit)
  expect_false(out$meta$lethal)
})

test_that("extract_tdt with lethal = TRUE returns a sensible T_crit", {
  skip_unless_brms()

  wf  <- load_fixture_workflow()
  ts  <- truth_summary()
  out <- suppressMessages(
    extract_tdt(wf, t_ref = 60, ndraws = 500, lethal = TRUE)
  )

  expect_true(!is.null(out$T_crit))
  expect_true(out$meta$lethal)

  # T_crit (rate-multiplier-based) sits below CTmax by 2-3 z; with default
  # TC_rate_range = c(0.1, 1) the median is at CTmax - 2.5 * z (geometric mean
  # of the log10(rate) range). The true T_crit median should fall inside the
  # posterior 95% CrI.
  true_T_crit_median <- ts$CTmax_1hr - 2.5 * ts$z
  expect_lt(out$T_crit$summary$temp_median, out$CTmax$summary$temp_median)
  expect_lte(out$T_crit$summary$temp_lower, true_T_crit_median)
  expect_gte(out$T_crit$summary$temp_upper, true_T_crit_median)

  # T_crit's 95% CrI should be approximately as wide as z (since pooling
  # uniformly over log10(r*) ∈ [-3, -2] gives a range of 1 z, with parameter
  # uncertainty adding a smaller component). Sanity check: width > 0.5 z.
  ci_width <- out$T_crit$summary$temp_upper - out$T_crit$summary$temp_lower
  expect_gt(ci_width, 0.5 * ts$z)
})

test_that("extract_tdt z is read from the posterior and z_local is opt-in", {
  skip_unless_brms()

  wf <- load_fixture_workflow()

  # Default: no local-z block.
  out <- extract_tdt(wf, t_ref = 60, ndraws = 500)
  expect_null(out$z$local)
  expect_named(out$z, c("draws", "summary", "local"), ignore.order = TRUE)

  # Relative z equals -1 / b_mid_temp_c read straight from the posterior — no
  # regression. (Compare medians over all draws; extract_tdt subsamples.)
  d  <- posterior::as_draws_df(get_brmsfit(wf)) |> as.data.frame()
  out_all <- extract_tdt(wf, t_ref = 60, target_surv = "relative",
                         ndraws = brms::ndraws(get_brmsfit(wf)))
  expect_equal(out_all$z$summary$z_median,
               stats::median(-1 / d$b_mid_temp_c), tolerance = 1e-8)

  # z_local = TRUE returns per-temperature local z(T).
  out_loc <- extract_tdt(wf, t_ref = 60, ndraws = 500,
                         target_surv = "absolute", z_local = TRUE)
  expect_false(is.null(out_loc$z$local))
  expect_true(all(c("temp", "z_median", "z_lower", "z_upper") %in%
                  names(out_loc$z$local$summary)))
  expect_true(all(is.finite(out_loc$z$local$summary$z_median)))
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

test_that("an absolute target_surv outside the fitted asymptotes warns clearly", {
  skip_unless_brms()
  wf <- load_fixture_workflow_grouped()   # survival spans ~[0.04, 0.97]
  # 0.1 lies below the lower asymptote for some draws -> LT undefined there:
  # a clear, reasoned warning instead of a bare "NaNs produced".
  expect_warning(extract_tdt(wf, target_surv = 0.1, ndraws = NULL),
                 "outside the fitted curve's asymptotes")
  # an in-range absolute threshold does not raise that warning
  asymp_warn <- FALSE
  withCallingHandlers(
    extract_tdt(wf, target_surv = 0.5, ndraws = NULL),
    warning = function(w) {
      if (grepl("asymptotes", conditionMessage(w))) asymp_warn <<- TRUE
      invokeRestart("muffleWarning")
    })
  expect_false(asymp_warn)
})

test_that("extract_tdt clamps ndraws to the posterior size (no crash when ndraws > draws)", {
  skip_unless_brms()
  wf <- load_fixture_workflow()                 # ~1500 posterior draws
  # Regression: relative-mode derive_tdt_curve() passed `ndraws` straight to
  # brms::posterior_linpred(), which errors when ndraws exceeds the posterior
  # size. Requesting more draws than exist must now clamp, not crash. This is
  # also what the default ndraws = 1000 hits on a small fit (< 1000 draws).
  expect_no_error(extract_tdt(wf, target_surv = "relative", ndraws = 5000))
  expect_no_error(derive_tdt_curve(wf, temp_grid = c(30, 34, 38), ndraws = 5000))
})
