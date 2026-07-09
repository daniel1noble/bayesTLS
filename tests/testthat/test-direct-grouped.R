# Gated per-group tests for the group-aware extraction engine. A 2-group direct
# fit (ctmax/z ~ 0 + grp) with a deliberately large per-group z gap (A z=5.56,
# B z=3.33). Confirms: (a) per-group truth recovery via tls(by=); (b)
# coding-invariance (~ 0 + grp vs ~ grp); (c) the by= predictors return
# per-group output. Set RUN_BRMS_TESTS=true.

test_that("tls(by=) recovers per-group truth with the correct group->z mapping", {
  skip_unless_brms()
  wf  <- load_fixture_workflow_grouped()
  tr  <- attr(wf, "truth_2group")
  out <- tls(wf, by = "grp", params = c("z", "ctmax"), t_ref = 60, ndraws = 1000,
             temp_mean = wf$meta$temp_mean)$summary
  gg  <- function(g, q) out$median[out$grp == g & out$quantity == q]

  expect_true("grp" %in% names(out))                       # per-group output
  expect_setequal(out$grp, c("A", "B"))
  zA <- gg("A", "z"); zB <- gg("B", "z")
  cA <- gg("A", "CTmax"); cB <- gg("B", "CTmax")
  # per-group z recovers truth within tolerance and the A>B ordering is correct
  expect_lt(abs(zA - tr$A$z), 1.5)
  expect_lt(abs(zB - tr$B$z), 1.5)
  expect_gt(zA - zB, 1.0)                                   # large gap -> mapping is real
  expect_lt(abs(cA - tr$A$CTmax_1hr), 1.0)
  expect_lt(abs(cB - tr$B$CTmax_1hr), 1.0)
})

test_that("per-group z is invariant to cell-means (~ 0 + grp) vs treatment (~ grp) coding", {
  skip_unless_brms()
  wf_cm <- load_fixture_workflow_grouped()                 # ~ 0 + grp
  raw <- simulate_tdt_2group(seed = 20260622)
  std <- standardize_data(raw, temp = "T", duration = "t_hours",
                          n_total = "n_trials", n_surv = "y_alive", duration_unit = "hours")
  wf_tr <- fit_4pl(std, ctmax = ~ grp, z = ~ grp,          # treatment coding
                   chains = 2, iter = 1200, warmup = 600, cores = 2, seed = 1,
                   control = list(adapt_delta = 0.95, max_treedepth = 12),
                   file = file.path(here::here("tests", "testthat", "fixtures"),
                                    "sim_fit_grouped_trt_small"), refresh = 0)
  z_cm <- tls(wf_cm, params = "z", by = "grp", ndraws = NULL,
              temp_mean = wf_cm$meta$temp_mean)$summary
  z_tr <- tls(wf_tr, params = "z", by = "grp", ndraws = NULL,
              temp_mean = wf_tr$meta$temp_mean)$summary
  gm <- function(s, g) s$median[s$grp == g & s$quantity == "z"]
  # The two codings share the LIKELIHOOD and, since make_4pl_priors() routes
  # CTmaxdev/logz through a coding-aware builder (contrasts centred on ZERO, not
  # log(3)), they also share the effective prior on every per-group z. Per-group z
  # must therefore agree up to MCMC noise across two fresh fits. Validated: |diff|
  # = 0.006 (A) / 0.144 (B) post-fix. Before the coding-aware prior fix this was
  # ~2.5 (group-B z 6.7 cell-means vs 9.2 treatment); tol = 0.3 catches any
  # regression of that bug.
  for (gg in c("A", "B")) {
    expect_lt(abs(gm(z_cm, gg) - gm(z_tr, gg)), 0.3)
  }
})

test_that("plot helpers facet grouped output instead of pooling across groups", {
  skip_unless_brms()
  wf <- load_fixture_workflow_grouped()
  sc <- predict_survival_curves(wf, temps = c(33, 36), durations = c(1, 8), ndraws = 100)
  p_sc <- plot_survival_curves(sc)
  expect_s3_class(p_sc, "ggplot")
  expect_s3_class(p_sc$facet, "FacetWrap")                 # per-group panels
  tg <- seq(min(wf$data$temp) - 2, max(wf$data$temp) + 2, by = 0.1)
  td <- plot_temperature_density(
    derive_temperature_for_duration(wf, exposure_duration = 1, temp_grid = tg,
                                    by = "grp", seed = 1))
  expect_s3_class(td$facet, "FacetWrap")
  ltx <- derive_tdt_curve(wf, temp_grid = c(33, 36), by = "grp", ndraws = 100)
  expect_true("grp" %in% names(ltx$summary))
  expect_setequal(unique(ltx$summary$grp), c("A", "B"))
  p_ltx <- plot_tdt_curve(ltx, panels = "log")
  expect_s3_class(p_ltx$facet, "FacetWrap")
})

test_that("grouped derived-quantity and survival draws carry grp (no .draw key collision)", {
  skip_unless_brms()
  wf <- load_fixture_workflow_grouped()
  dd <- tls(wf, by = "grp", lethal = TRUE, ndraws = 300, seed = 1,
            temp_mean = wf$meta$temp_mean)$draws
  expect_true("grp" %in% names(dd)); expect_setequal(unique(dd$grp), c("A", "B"))
  # (grp, quantity, .draw) keys are unique -> groups not collapsed onto colliding .draw
  expect_equal(nrow(dplyr::distinct(dd, grp, quantity, .draw)), nrow(dd))

  sd <- get_surv_draws(predict_survival_curves(wf, temps = c(33, 36),
                                               durations = c(1, 8), ndraws = 100))
  expect_true("grp" %in% names(sd))
  expect_equal(nrow(dplyr::distinct(sd, grp, temp, duration, .draw)), nrow(sd))

  ltx <- derive_tdt_curve(wf, temp_grid = c(33, 36), by = "grp", ndraws = 100)
  expect_true("grp" %in% names(ltx$draws))
  expect_equal(nrow(dplyr::distinct(ltx$draws, grp, temp, target_surv, .draw)),
               nrow(ltx$draws))
})

test_that("predict_survival_curves(by=) and predict_heat_injury(by=) return per-group output", {
  skip_unless_brms()
  wf <- load_fixture_workflow_grouped()
  sc <- predict_survival_curves(wf, temps = c(33, 36), durations = c(1, 8), ndraws = 300)$summary
  expect_true("grp" %in% names(sc))
  expect_setequal(unique(sc$grp), c("A", "B"))
  expect_true(all(is.finite(sc$survival_median)))

  hi <- suppressMessages(predict_heat_injury(
    data.frame(time = 0:5, temp = seq(33, 38, length.out = 6)), wf, ndraws = 300, seed = 1))$summary
  expect_true("grp" %in% names(hi))
  expect_setequal(unique(hi$grp), c("A", "B"))
  expect_true(all(is.finite(hi$surv_median)))
})
