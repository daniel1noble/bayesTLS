test_that("make_temperature_scenarios returns four traces of the right shape", {
  s <- make_temperature_scenarios(n_hours = 48, baseline = 20, spike_temp = 30,
                                  diurnal_n_days = 2)
  expect_named(s, c("flat", "single_spike", "multi_spike", "diurnal"))
  for (nm in c("flat", "single_spike", "multi_spike")) {
    expect_s3_class(s[[nm]], "tbl_df")
    expect_equal(nrow(s[[nm]]), 48)
    expect_named(s[[nm]], c("time", "temp"))
  }
  expect_s3_class(s$diurnal, "tbl_df")
  expect_equal(nrow(s$diurnal), 2 * 24)
  expect_named(s$diurnal, c("time", "temp"))
})

test_that("diurnal night-temp vector lets cool/warm days have distinct night baselines", {
  s <- make_temperature_scenarios(
    n_hours = 96, baseline = 18, spike_temp = 30,
    diurnal_n_days     = 2,
    diurnal_night_temp = c(15, 22),    # cool night day 1, warm night day 2
    diurnal_day_peaks  = c(20, 30),
    diurnal_noise_sd   = 0
  )
  d <- s$diurnal
  # 02:00 sits in the cool tail of each day, far from the 14:00 Gaussian peak,
  # so the trace temp there should be very close to that day's night baseline.
  night_d1 <- d$temp[d$time == 2]
  night_d2 <- d$temp[d$time == 26]
  expect_equal(night_d1, 15, tolerance = 0.05)
  expect_equal(night_d2, 22, tolerance = 0.05)
  # Sanity: hot-day peak hotter than cool-day peak
  peak_d1 <- max(d$temp[d$time < 24])
  peak_d2 <- max(d$temp[d$time >= 24])
  expect_gt(peak_d2, peak_d1)
})

test_that("flat trace is constant baseline; spike traces inject spikes at the right hours", {
  s <- make_temperature_scenarios(baseline = 20, spike_temp = 30,
                                  n_hours = 48,
                                  spike_times_single = 12,
                                  spike_times_multi  = c(12, 24, 36))
  expect_true(all(s$flat$temp == 20))
  expect_equal(sum(s$single_spike$temp == 30), 1)
  expect_equal(s$single_spike$temp[s$single_spike$time == 12], 30)
  expect_equal(sum(s$multi_spike$temp == 30), 3)
})

test_that("planted_dose_from_trace returns zero for a sub-T_c flat trace", {
  s <- make_temperature_scenarios(baseline = 20, spike_temp = 30,
                                  n_hours = 48)
  d <- planted_dose_from_trace(s$flat, z = 5, CTmax_1hr = 32, T_c = 25)
  expect_equal(d$hi_inc, rep(0, nrow(s$flat)))
  expect_equal(tail(d$hi_cumulative, 1), 0)
})

test_that("planted_dose_from_trace matches the analytical formula for a single spike", {
  z   <- 5
  CT1 <- 32
  T_c <- 25
  s_temp <- 30
  s   <- make_temperature_scenarios(baseline = 20, spike_temp = s_temp,
                                    n_hours = 48,
                                    spike_times_single = 12)
  d   <- planted_dose_from_trace(s$single_spike, z = z, CTmax_1hr = CT1, T_c = T_c)
  expected <- 100 * 10 ^ ((s_temp - CT1) / z)
  expect_equal(tail(d$hi_cumulative, 1), expected, tolerance = 1e-8)
})

test_that("single-spike dose scales with dt_hours (one-time-step spike width)", {
  z   <- 5
  CT1 <- 32
  T_c <- 25
  s_temp <- 30
  dt  <- 0.5
  s   <- make_temperature_scenarios(baseline = 20, spike_temp = s_temp,
                                     n_hours = 48, dt_hours = dt,
                                     spike_times_single = 12)
  d   <- planted_dose_from_trace(s$single_spike, z = z, CTmax_1hr = CT1, T_c = T_c)
  # Each spike is one time step (dt_hours) wide, so the delivered dose is the
  # instantaneous rate times dt_hours (documented behaviour).
  expected <- 100 * 10 ^ ((s_temp - CT1) / z) * dt
  expect_equal(tail(d$hi_cumulative, 1), expected, tolerance = 1e-8)
})

test_that("planted_dose_from_trace cumulative equals sum of single-spike contributions", {
  z   <- 5
  CT1 <- 32
  T_c <- 25
  s_temp <- 30
  spike_hours <- c(12, 24, 36)
  s   <- make_temperature_scenarios(baseline = 20, spike_temp = s_temp,
                                    n_hours = 48,
                                    spike_times_multi = spike_hours)
  d   <- planted_dose_from_trace(s$multi_spike, z = z, CTmax_1hr = CT1, T_c = T_c)
  expected_total <- length(spike_hours) * 100 * 10 ^ ((s_temp - CT1) / z)
  expect_equal(tail(d$hi_cumulative, 1), expected_total, tolerance = 1e-8)
})
