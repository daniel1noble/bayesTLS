# `local` gate in tdt_z_from_pars(): skipping the per-temperature breakdown must
# not change the pooled z, and must only suppress local_draws/local_summary.
# Uses synthetic coefficient draws so no brms fit is required (fast unit test).

make_pars <- function(np = 300, seed = 42) {
  set.seed(seed)
  data.frame(
    b_mid_Intercept    = stats::rnorm(np,  1.5,  0.05),
    b_mid_temp_c       = stats::rnorm(np, -0.15, 0.01),  # z = -1/slope > 0
    b_lowraw_Intercept = stats::rnorm(np, -2.0,  0.10),  # low well below 0.5
    b_lowraw_temp_c    = stats::rnorm(np,  0.00, 0.02),
    b_upraw_Intercept  = stats::rnorm(np,  1.5,  0.10),  # up well above 0.5
    b_upraw_temp_c     = stats::rnorm(np, -0.05, 0.01),  # up declines with T -> bends LT
    b_logk_Intercept   = stats::rnorm(np, log(8), 0.10),
    b_logk_temp_c      = stats::rnorm(np,  0.02, 0.01)
  )
}

bnd  <- list(low_min = 0.001, low_w = 0.498, up_min = 0.501, up_w = 0.498)
Tbar <- 34.5
grid <- c(30, 33, 36, 39)
ts_rel <- list(mode = "relative", prob = NA_real_, label = "(low+up)/2")
ts_abs <- list(mode = "absolute", prob = 0.5,      label = "p=0.500")

for (nm in c("relative", "absolute")) {
  ts <- if (nm == "relative") ts_rel else ts_abs

  test_that(sprintf("[%s] pooled z is identical with local on vs off", nm), {
    pars <- make_pars()
    z_off <- tdt_z_from_pars(pars, Tbar, bnd, ts, grid, local = FALSE)
    z_on  <- tdt_z_from_pars(pars, Tbar, bnd, ts, grid, local = TRUE)
    # The reported z (draws + summary) must be byte-identical regardless of the gate.
    expect_identical(z_off$draws,   z_on$draws)
    expect_identical(z_off$summary, z_on$summary)
  })

  test_that(sprintf("[%s] local block is NULL only when local = FALSE", nm), {
    pars <- make_pars()
    z_off <- tdt_z_from_pars(pars, Tbar, bnd, ts, grid, local = FALSE)
    z_on  <- tdt_z_from_pars(pars, Tbar, bnd, ts, grid, local = TRUE)
    expect_null(z_off$local_draws)
    expect_null(z_off$local_summary)
    expect_s3_class(z_on$local_draws,   "tbl_df")
    expect_s3_class(z_on$local_summary, "tbl_df")
    # One summary row per finite temperature on the grid.
    expect_equal(nrow(z_on$local_summary), length(grid))
  })

  test_that(sprintf("[%s] default is local = TRUE (preserves derive_z contract)", nm), {
    pars <- make_pars()
    z_def <- tdt_z_from_pars(pars, Tbar, bnd, ts, grid)  # no local arg
    expect_s3_class(z_def$local_summary, "tbl_df")
  })
}

test_that("absolute local z genuinely varies across temperature (gate exercises the bent-curve path)", {
  # Sanity: with up & k temperature-dependent, the per-temperature local z is not
  # constant, so the local block carries real (skippable) information.
  pars <- make_pars()
  z_on <- tdt_z_from_pars(pars, Tbar, bnd, ts_abs, grid, local = TRUE)
  expect_gt(stats::sd(z_on$local_summary$z_median), 0)
})
