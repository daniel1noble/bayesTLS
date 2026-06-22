# Gated per-group tests for the group-aware extraction engine. A 2-group direct
# fit (ctmax/z ~ 0 + grp) with a deliberately large per-group z gap (A z=5.56,
# B z=3.33). Confirms: (a) per-group truth recovery; (b) extract_tdt(by=) agrees
# with tls(by=) per group; (c) coding-invariance (~ 0 + grp vs ~ grp); (d) the
# by= predictors return per-group output. Set RUN_BRMS_TESTS=true.

test_that("extract_tdt(by=) recovers per-group truth with the correct group->z mapping", {
  skip_unless_brms()
  wf  <- load_fixture_workflow_grouped()
  tr  <- attr(wf, "truth_2group")
  out <- extract_tdt(wf, t_ref = 60, ndraws = 1000)        # by = NULL -> auto group by grp

  expect_true("grp" %in% names(out$z$summary))             # per-group output
  expect_setequal(out$z$summary$grp, c("A", "B"))
  zA <- out$z$summary$z_median[out$z$summary$grp == "A"]
  zB <- out$z$summary$z_median[out$z$summary$grp == "B"]
  cA <- out$CTmax$summary$temp_median[out$CTmax$summary$grp == "A"]
  cB <- out$CTmax$summary$temp_median[out$CTmax$summary$grp == "B"]
  # per-group z recovers truth within tolerance and the A>B ordering is correct
  expect_lt(abs(zA - tr$A$z), 1.5)
  expect_lt(abs(zB - tr$B$z), 1.5)
  expect_gt(zA - zB, 1.0)                                   # large gap -> mapping is real
  expect_lt(abs(cA - tr$A$CTmax_1hr), 1.0)
  expect_lt(abs(cB - tr$B$CTmax_1hr), 1.0)
})

test_that("extract_tdt(by=) agrees with tls(by=) per group", {
  skip_unless_brms()
  wf  <- load_fixture_workflow_grouped()
  et  <- extract_tdt(wf, t_ref = 60, ndraws = NULL)
  tl  <- tls(wf, by = "grp", params = c("z", "ctmax"), mode = "relative",
             t_ref = 60, temp_mean = wf$meta$temp_mean)$summary
  for (gg in c("A", "B")) {
    z_et <- et$z$summary$z_median[et$z$summary$grp == gg]
    z_tl <- tl$median[tl$grp == gg & tl$quantity == "z"]
    expect_lt(abs(z_et - z_tl), 0.3)                        # two engines, same model
  }
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
  z_cm <- extract_tdt(wf_cm, ndraws = NULL)$z$summary
  z_tr <- extract_tdt(wf_tr, ndraws = NULL)$z$summary
  # The two codings share the LIKELIHOOD but not the prior parameterisation
  # (cell-means: one prior per level; treatment: prior on the reference + on each
  # contrast, so group B gets a wider prior), so on small fits per-group z differs
  # by a few tenths from prior parameterisation + MCMC noise across two fresh fits.
  # The real failure this guards is the treatment-coding reference-level bug, which
  # would report group B's z AS group A's -> a |z_A - z_B| ~ 2.2 gap; tol = 1.0
  # tolerates the benign difference while still catching that.
  for (gg in c("A", "B")) {
    expect_lt(abs(z_cm$z_median[z_cm$grp == gg] - z_tr$z_median[z_tr$grp == gg]), 1.0)
  }
})

test_that("plot helpers facet grouped output instead of pooling across groups", {
  skip_unless_brms()
  wf <- load_fixture_workflow_grouped()
  sc <- predict_survival_curves(wf, temps = c(33, 36), durations = c(1, 8), ndraws = 100)
  p_sc <- plot_survival_curves(sc)
  expect_s3_class(p_sc, "ggplot")
  expect_s3_class(p_sc$facet, "FacetWrap")                 # per-group panels
  td <- plot_temperature_density(extract_tdt(wf, ndraws = 200, seed = 1)$CTmax)
  expect_s3_class(td$facet, "FacetWrap")
})

test_that("accessors carry the group column on grouped output (no .draw key collision)", {
  skip_unless_brms()
  wf <- load_fixture_workflow_grouped()
  et <- extract_tdt(wf, ndraws = 300, lethal = TRUE, seed = 1)
  zc <- get_ctmax_draws(et)
  tc <- get_tcrit_draws(et)
  expect_true("grp" %in% names(zc)); expect_setequal(unique(zc$grp), c("A", "B"))
  expect_true("grp" %in% names(tc)); expect_setequal(unique(tc$grp), c("A", "B"))
  # (grp, .draw) keys are unique -> groups not collapsed onto colliding .draw
  expect_equal(nrow(dplyr::distinct(zc, grp, .draw)), nrow(zc))

  sd <- get_surv_draws(predict_survival_curves(wf, temps = c(33, 36),
                                               durations = c(1, 8), ndraws = 100))
  expect_true("grp" %in% names(sd))
  expect_equal(nrow(dplyr::distinct(sd, grp, temp, duration, .draw)), nrow(sd))
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
