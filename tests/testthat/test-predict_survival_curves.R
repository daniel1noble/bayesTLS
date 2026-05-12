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
