test_that("threshold_x_by_draw recovers the exact crossing on monotone toy data", {
  # One draw with a perfect 4PL-shaped row: survival decreasing in duration.
  # The crossing at survival = 0.5 should fall at duration = 2.
  x  <- c(0.5, 1, 2, 4, 8)
  ps <- c(0.99, 0.90, 0.50, 0.10, 0.01)
  pm <- matrix(ps, nrow = 1)
  expect_equal(threshold_x_by_draw(pm, x, target = 0.5), 2)
})

test_that("threshold_x_by_draw returns NA when target outside the range", {
  x  <- c(0.5, 1, 2, 4, 8)
  ps <- c(0.99, 0.90, 0.50, 0.10, 0.01)
  pm <- matrix(ps, nrow = 1)
  expect_true(is.na(threshold_x_by_draw(pm, x, target = 1.05)))
  expect_true(is.na(threshold_x_by_draw(pm, x, target = -0.01)))
})

test_that("threshold_x_by_draw handles a vector of draws", {
  # Two draws, both monotone decreasing but with different crossings.
  x   <- c(0.5, 1, 2, 4, 8)
  pm  <- rbind(
    c(0.99, 0.90, 0.50, 0.10, 0.01),  # crosses 0.5 at x = 2
    c(0.99, 0.95, 0.75, 0.50, 0.10)   # crosses 0.5 at x = 4
  )
  res <- threshold_x_by_draw(pm, x, target = 0.5)
  expect_length(res, 2)
  expect_equal(res, c(2, 4))
})

test_that("threshold_x_by_draw returns NA for a draw with non-finite predictions", {
  x  <- c(0.5, 1, 2, 4, 8)
  pm <- rbind(c(0.99, NA, 0.50, 0.10, 0.01),
              c(0.99, 0.90, 0.50, 0.10, 0.01))
  res <- threshold_x_by_draw(pm, x, target = 0.5)
  expect_true(is.na(res[1]))
  expect_equal(res[2], 2)
})

# --- summarise_observed_survival: pure data summary, no fit needed ----------

test_that("summarise_observed_survival computes per-cell mean, SE and clamped bounds", {
  # Two cells with known survival values; verify against hand computation.
  obs <- tibble::tibble(
    temp     = c(30, 30, 30, 34, 34, 34),
    duration = c(1, 1, 1, 2, 2, 2),
    survival = c(0.9, 0.8, 1.0, 0.2, 0.1, 0.0),
    n_total  = c(20, 20, 20, 20, 20, 20))
  out <- summarise_observed_survival(obs)

  expect_equal(nrow(out), 2L)
  expect_true(all(c("survival_mean", "survival_se", "survival_lower",
                    "survival_upper", "n_units", "n_total_sum") %in% names(out)))

  cell30 <- out[out$temp == 30, ]
  v30 <- c(0.9, 0.8, 1.0)
  expect_equal(cell30$survival_mean, mean(v30))
  expect_equal(cell30$survival_se,   stats::sd(v30) / sqrt(3))
  expect_equal(cell30$n_units, 3L)
  expect_equal(cell30$n_total_sum, 60)

  # Upper bound for the 30 C cell would exceed 1 (mean 0.9 + SE) and is clamped.
  expect_lte(cell30$survival_upper, 1)
  expect_equal(cell30$survival_upper, min(1, cell30$survival_mean + cell30$survival_se))
  # Lower bound for the 34 C cell would dip below 0 and is clamped at 0.
  cell34 <- out[out$temp == 34, ]
  expect_gte(cell34$survival_lower, 0)
})

test_that("summarise_observed_survival drops n_total_sum for continuous-proportion data", {
  obs <- tibble::tibble(temp = c(30, 30), duration = c(1, 1),
                        survival = c(0.7, 0.5))     # no n_total column (Beta)
  out <- summarise_observed_survival(obs)
  expect_false("n_total_sum" %in% names(out))
  expect_equal(out$survival_mean, 0.6)
})

# --- predict_survival_curves: gated, against the fixture's known truth ------

test_that("predict_survival_curves returns a valid, monotone surface with asymptotes near truth", {
  skip_unless_brms()
  wf <- load_fixture_workflow()

  # Default temps are the unique assay temperatures.
  pred <- predict_survival_curves(wf, durations = c(1, 5, 50), ndraws = 300)
  expect_setequal(unique(pred$summary$temp), sort(unique(wf$data$temp)))
  expect_equal(nrow(pred$summary), length(unique(wf$data$temp)) * 3L)
  expect_equal(dim(pred$draws_matrix), c(300L, nrow(pred$summary)))

  s <- pred$summary
  expect_true(all(s$survival_median >= 0 & s$survival_median <= 1))
  expect_true(all(s$survival_lower <= s$survival_median + 1e-9 &
                  s$survival_median <= s$survival_upper + 1e-9))

  # At a fixed temperature survival falls with duration.
  t1 <- s[s$temp == sort(unique(s$temp))[1], ]
  t1 <- t1[order(t1$duration), ]
  expect_true(all(diff(t1$survival_median) <= 1e-8))

  # Asymptote behaviour: very short exposure -> near the upper asymptote
  # (truth u = 0.95), very long exposure -> near the lower asymptote
  # (truth ell = 0.05). Evaluate at the coolest assay temperature.
  t_cool <- sort(unique(wf$data$temp))[1]
  ends   <- predict_survival_curves(wf, temps = t_cool,
                                    durations = c(1e-3, 1e4), ndraws = 300)$summary
  expect_gt(ends$survival_median[ends$duration == min(ends$duration)], 0.8)
  expect_lt(ends$survival_median[ends$duration == max(ends$duration)], 0.2)
})
