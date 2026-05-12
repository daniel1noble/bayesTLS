test_that("time_to_surv_threshold_4pl returns NA when target outside (low, up)", {
  # target = 0.05 is below low = 0.1, so the 4PL never reaches it → NA
  expect_true(is.na(time_to_surv_threshold_4pl(
    temp = 30, survival_target = 0.05,
    low = 0.1, up = 0.95, k = 6,
    mid_int = 1, mid_temp = -0.2, temp_mean = 30
  )))
})

test_that("time_to_surv_threshold_4pl recovers the analytical inverse at T_bar", {
  # At T_bar (temp_c = 0), mid = mid_int. Survival = 0.5 at duration = 10^mid_int
  # when (low, up) are symmetric around 0.5.
  out <- time_to_surv_threshold_4pl(
    temp = 30, survival_target = 0.5,
    low = 0.05, up = 0.95, k = 6,
    mid_int = 1.5, mid_temp = -0.2, temp_mean = 30
  )
  # mid = 1.5, log term = log(0.45/0.45)/k = 0, so duration = 10^1.5
  expect_equal(out, 10 ^ 1.5, tolerance = 1e-6)
})

test_that("survival_from_dose returns u at dose = 0 and target_surv at dose = 1", {
  expect_equal(survival_from_dose(0,
                                  low = 0.02, up = 0.98, k = 6,
                                  target_surv = 0.5),
               0.98, tolerance = 1e-6)
  expect_equal(survival_from_dose(1,
                                  low = 0.02, up = 0.98, k = 6,
                                  target_surv = 0.5),
               0.5, tolerance = 1e-6)
})

test_that("repair_rate_schoolfield returns positive rates peaked near TREF", {
  rp <- list(TA = 14065, TAL = 50000, TAH = 120000,
             TL = 10.5 + 273.15, TH = 22.5 + 273.15,
             TREF = 17 + 273.15, r_ref = 0.01)
  temps <- seq(5, 30, by = 1)
  rates <- repair_rate_schoolfield(
    temp_celsius = temps,
    TA = rp$TA, TAL = rp$TAL, TAH = rp$TAH,
    TL = rp$TL, TH = rp$TH, TREF = rp$TREF,
    r_ref = rp$r_ref
  )
  expect_true(all(rates >= 0))
  expect_equal(length(rates), length(temps))
  # The peak should be near TREF = 17 °C (TPC has its optimum in this region).
  peak <- temps[which.max(rates)]
  expect_true(peak >= 13 && peak <= 22)
})
