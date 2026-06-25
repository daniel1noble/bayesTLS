# Tests for derive_tdt_landscape() (R/tdt_landscape.R). The function builds a
# dense default temperature x duration grid and delegates to
# predict_survival_curves(). We check the default grid construction, the
# physical constraints on the returned surface, and that it delegates without
# altering the grid (identical output to a direct predict_survival_curves()
# call on the same grid). Gated behind RUN_BRMS_TESTS.

test_that("derive_tdt_landscape rejects an unfitted workflow", {
  expect_error(derive_tdt_landscape(list(fit = NULL)), "NULL")
})

test_that("default duration grid is linearly spaced (regular for geom_raster)", {
  # Fast-tier regression for the uneven-raster bug: the default duration grid
  # must be evenly spaced on the LINEAR axis, not log-spaced. Stub out the fit
  # check and the delegate so no brms model is needed; just capture the grid
  # derive_tdt_landscape() constructs.
  wf <- list(fit = TRUE, data = data.frame(temp = c(30, 40),
                                            duration = c(1, 60)))
  captured <- NULL
  testthat::with_mocked_bindings(
    has_fit = function(...) TRUE,
    predict_survival_curves = function(workflow, temps, durations, ...) {
      captured <<- durations
      list(summary = NULL)
    },
    {
      derive_tdt_landscape(wf)
    }
  )
  expect_length(captured, 120L)
  expect_equal(range(captured), c(1, 60), tolerance = 1e-8)
  expect_equal(stats::sd(diff(captured)), 0, tolerance = 1e-8)   # even (linear)
})

test_that("derive_tdt_landscape builds a 120x120 grid spanning the data range", {
  skip_unless_brms()
  wf  <- load_fixture_workflow()
  lsp <- derive_tdt_landscape(wf, ndraws = 200)

  expect_true(all(c("temp", "duration", "survival_median",
                    "survival_lower", "survival_upper") %in% names(lsp$summary)))
  expect_equal(nrow(lsp$summary), 120L * 120L)         # default 120 x 120

  # Default temperature grid spans the observed temperature range exactly.
  trange <- range(wf$data$temp)
  expect_equal(range(lsp$summary$temp), trange, tolerance = 1e-8)
  # Default duration grid spans the observed duration range and is evenly
  # spaced on the LINEAR axis (regular grid for geom_raster in
  # plot_tdt_landscape, which uses a linear duration axis by default).
  drange <- range(wf$data$duration)
  expect_equal(range(lsp$summary$duration), drange, tolerance = 1e-6)
  ud <- sort(unique(lsp$summary$duration))
  expect_equal(stats::sd(diff(ud)), 0, tolerance = 1e-6)  # even (linear)
})

test_that("derive_tdt_landscape returns valid, monotone survival probabilities", {
  skip_unless_brms()
  wf  <- load_fixture_workflow()
  lsp <- derive_tdt_landscape(wf, ndraws = 200)

  s <- lsp$summary
  expect_true(all(s$survival_median >= 0 & s$survival_median <= 1))
  expect_true(all(s$survival_lower <= s$survival_median + 1e-9))
  expect_true(all(s$survival_median <= s$survival_upper + 1e-9))

  # Survival falls with exposure duration at a fixed temperature.
  one_temp <- s[s$temp == sort(unique(s$temp))[6], ]
  one_temp <- one_temp[order(one_temp$duration), ]
  expect_true(all(diff(one_temp$survival_median) <= 1e-8))
})

test_that("derive_tdt_landscape delegates to predict_survival_curves unchanged", {
  skip_unless_brms()
  wf <- load_fixture_workflow()
  tg <- c(31, 34, 37)
  dg <- c(0.5, 5, 50)
  # Using ALL posterior draws makes both calls deterministic, so the surfaces
  # must match exactly (the landscape is predict_survival_curves on this grid).
  nd  <- posterior::ndraws(posterior::as_draws(wf$fit))
  lsp <- derive_tdt_landscape(wf, temp_grid = tg, duration_grid = dg, ndraws = nd)
  psc <- predict_survival_curves(wf, temps = tg, durations = dg, ndraws = nd)
  expect_equal(lsp$summary, psc$summary)
})
