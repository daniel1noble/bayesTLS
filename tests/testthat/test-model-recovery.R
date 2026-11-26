# Alternative-models Phase 3: end-to-end recovery. Data are simulated under
# each dose-response model with known low, up, k, z and CTmax
# (simulate_tdt(model =)) in the Phase 3 design (helper-skip.R: 30, 33, 36 degC
# x 7 durations every half decade from 0.05 to 50 h x 3 replicates of 30,
# beta-binomial noise) and fitted with the matching model (cached fixtures,
# load_fixture_workflow_model()). The design matters. With durations a decade
# apart a steep curve can slide each temperature's midpoint between two
# durations, the data barely pin z, and recovery checks pass whatever the code
# does to z (exact large-sample comparison of designs in
# notes/package_design/alt_models_phase3/R/kl_design_check.R).
#
# Two kinds of check. Recovery of the truth: point tolerances that cover every
# matching fit in ten seeds of this design (alt_models_phase3/R/
# phase3_seed_check.R), and the truth inside a wide 99.8% interval (a
# calibrated 95% interval misses its truth 5% of the time, so across dozens of
# such checks some would fail by chance; with 1,200 draws the 99.8% bounds are
# set by the two or three most extreme draws, so they are approximate). And
# truth-free checks that tls() and tdt_parameter_table() report what the
# fitted model implies, draw by draw, which catch defects far smaller than any
# recovery tolerance. Fit-free versions of the per-draw checks are in
# test-models.R (swapped-formula Beta fixture; run by devtools::test(), not by
# R CMD check, whose tarball leaves the fixtures out). The fitted tests are
# gated behind RUN_BRMS_TESTS.

recovery_models <- c("log-logistic", "log-normal", "weibull1", "weibull2")

# Is `truth` inside the central 99.8% of the finite draws?
in_wide <- function(x, truth) {
  q <- stats::quantile(x[is.finite(x)], c(0.001, 0.999), names = FALSE)
  truth >= q[1] && truth <= q[2]
}

test_that("simulate_tdt() draws survival from the requested curve", {
  # each model's curve in its midpoint form, written out independently of the
  # registry; the truth is low 0.05, up 0.95, k 6
  G <- list("log-logistic" = function(u) stats::plogis(-u),
            "log-normal"   = function(u) stats::pnorm(-u),
            weibull1       = function(u) 2^(-exp(u)),
            weibull2       = function(u) 1 - 2^(-exp(-u)))
  ll <- simulate_tdt(seed = 7)
  for (m in names(G)) {
    d <- simulate_tdt(seed = 7, model = m)
    expect_identical(attr(d, "model"), m)
    expect_equal(d$p_true, 0.05 + 0.9 * G[[m]](6 * (d$log10_t_min - d$mid_t)),
                 tolerance = 1e-12, info = m)
    if (m != "log-logistic") expect_gt(max(abs(d$p_true - ll$p_true)), 0.03)
  }
})

test_that("each model recovers z and CTmax, relative and absolute, from data simulated under it", {
  skip_unless_brms()
  # Ten seeds of this design (phase3_seed_check.R), matching fits: the largest
  # errors were 1.8 in relative z, 2.5 in absolute-target z (absolute z is the
  # average slope of a fitted, slightly bent threshold line, so it varies
  # more), 0.6 degC in relative and 0.8 degC in absolute CTmax; every 99.8%
  # interval for z, CTmax and the gap below held its truth.
  z_tol <- c(relative = 2, `absolute 0.3` = 3)
  for (m in recovery_models) {
    wf <- load_fixture_workflow_model(m)
    tr <- attr(wf, "truth")
    expect_equal(wf$meta$model, m)
    # The absolute target 0.3 is read through each model's own correction, so
    # its true CTmax differs between models (relative CTmax does not).
    ct_abs <- truth_summary(tr$truth, model = m, p_abs = 0.3)$CTmax_1hr
    targets <- list(relative = list(ts = "relative", ct = tr$CTmax_rel_1hr),
                    `absolute 0.3` = list(ts = 0.3, ct = ct_abs))
    ct_draws <- list()
    for (nm in names(targets)) {
      tg  <- targets[[nm]]
      res <- tls(wf, params = c("z", "ctmax"), t_ref = 60, target_surv = tg$ts)
      s   <- res$summary
      dd  <- res$draws
      lab <- paste(m, nm)
      expect_lt(abs(s$median[s$quantity == "z"] - tr$z), z_tol[[nm]],
                label = paste(lab, "z"))
      expect_lt(abs(s$median[s$quantity == "CTmax"] - tg$ct), 1,
                label = paste(lab, "CTmax"))
      expect_true(in_wide(dd$value[dd$quantity == "z"], tr$z), info = lab)
      expect_true(in_wide(dd$value[dd$quantity == "CTmax"], tg$ct), info = lab)
      ct_draws[[nm]] <- dd[dd$quantity == "CTmax", c(".draw", "value")]
    }
    # The absolute threshold sits above the relative one by an amount set by
    # the curve and its asymptotes (0.55-0.88 degC here). Per draw, the gap
    # must bracket the true gap: this fails a reading that collapses the
    # absolute target onto the relative one. It is too wide to tell one
    # model's correction from another's; the per-draw test below does that.
    gap <- merge(ct_draws[["absolute 0.3"]], ct_draws[["relative"]], by = ".draw")
    expect_true(in_wide(gap$value.x - gap$value.y, ct_abs - tr$CTmax_rel_1hr), info = m)
  }
})

test_that("tls() CTmax and z match each fitted model's own threshold line, per draw", {
  skip_unless_brms()
  # References that never use the models' corrections or tls()'s handling of
  # t_ref. Relative target: the midpoint line brms evaluates
  # (derive_tdt_curve()), solved in temperature; this checks tls()'s crossing,
  # z and centring, exactly. Absolute target: brms's predicted survival solved
  # for the threshold time on a dense duration grid, then in temperature; its
  # z is the local slope of that line averaged over tls()'s default grid (11
  # temperatures across the data), as tls() averages it. CTmax is read at 60
  # and 30 min, so t_ref is not only ever the fit's own. These fits are in
  # hours.
  temps <- seq(30, 36, by = 0.25)
  ids   <- 1:300
  zg    <- seq(30, 36, length.out = 11)
  d     <- 0.1
  t_abs <- sort(unique(c(temps, zg - d, zg + d)))
  dgrid <- 10^seq(log10(0.02), log10(20), length.out = 400)
  pick  <- function(dd, q, draws) {
    x <- dd[dd$quantity == q, ]
    x$value[match(draws, x$.draw)]
  }
  for (m in recovery_models) {
    wf  <- load_fixture_workflow_model(m)
    fit <- get_brmsfit(wf)
    expect_identical(wf$meta$duration_unit, "hours")
    mid <- derive_tdt_curve(wf, temp_grid = temps)$draws
    M   <- matrix(NA_real_, max(mid$.draw), length(temps))
    M[cbind(mid$.draw, match(mid$temp, temps))] <- log10(mid$duration_model)
    g  <- new_tdt_grid(wf, temps = t_abs, durations = dgrid)
    pp <- brms::posterior_linpred(fit, newdata = g, transform = TRUE,
                                  re_formula = NA, draw_ids = ids)
    A  <- vapply(t_abs, function(t) {
      log10(threshold_x_by_draw(pp[, g$temp == t, drop = FALSE], dgrid, 0.3))
    }, numeric(length(ids)))
    on_t <- match(temps, t_abs)
    for (tr in c(60, 30)) {
      lab <- paste(m, "t_ref", tr)
      rel <- tls(wf, params = c("z", "ctmax"), t_ref = tr)$draws
      want <- tls_invert_logLT(M, log10(tr / 60), temps)
      got  <- pick(rel, "CTmax", seq_len(nrow(M)))
      ok   <- is.finite(want) & is.finite(got)
      expect_gt(mean(ok), 0.95)
      expect_equal(got[ok], want[ok], tolerance = 1e-8, info = lab)
      if (tr == 60) {
        expect_equal(pick(rel, "z", seq_len(nrow(M))),
                     -(temps[2] - temps[1]) / (M[, 2] - M[, 1]),
                     tolerance = 1e-8, info = lab)
      }
      ab   <- muffle_uncrossed(tls(wf, params = c("z", "ctmax"), t_ref = tr,
                                   target_surv = 0.3))$draws
      want <- tls_invert_logLT(A[, on_t], log10(tr / 60), temps)
      got  <- pick(ab, "CTmax", ids)
      ok   <- is.finite(want) & is.finite(got)
      expect_gt(mean(ok), 0.95)
      expect_lt(max(abs(got - want)[ok]), 0.02, label = paste(lab, "absolute CTmax"))
      if (tr == 60) {
        # a draw counts only if its threshold time is on the grid at all 22
        # temperatures (tls() averages all 11 local slopes)
        slope <- (A[, match(zg + d, t_abs)] - A[, match(zg - d, t_abs)]) / (2 * d)
        zref  <- rowMeans(-1 / slope)
        za    <- pick(ab, "z", ids)
        ok    <- is.finite(zref) & is.finite(za)
        expect_gt(mean(ok), 0.9)
        expect_lt(max(abs(za - zref)[ok]), 0.1, label = paste(m, "absolute z"))
      }
    }
  }
})

test_that("each model recovers its asymptotes and steepness at the centring temperature", {
  skip_unless_brms()
  for (m in recovery_models) {
    wf  <- load_fixture_workflow_model(m)
    tr  <- attr(wf, "truth")
    tab <- tdt_parameter_table(wf)
    expect_equal(attr(tab, "model"), m)
    # brms's draws of each sub-parameter at the centring temperature, put on
    # the natural scale by hand
    fit <- get_brmsfit(wf)
    nd  <- new_tdt_grid(wf, temps = wf$meta$temp_mean, durations = 1)
    lp  <- function(nlp) brms::posterior_linpred(fit, newdata = nd, nlpar = nlp,
                                                 re_formula = NA)[, 1]
    b   <- wf$meta$bounds
    draws <- list(low = b$low_min + stats::plogis(lp("lowraw")) * b$low_w,
                  up  = b$up_min + stats::plogis(lp("upraw")) * b$up_w,
                  k   = exp(lp("logk")))
    # The same asymptotes from brms's own predicted survival at a very short
    # and a very long exposure, independently of meta$bounds: so the bounds
    # written into the formula are the ones the workflow records.
    ends <- brms::posterior_linpred(fit, newdata = new_tdt_grid(wf, temps = wf$meta$temp_mean,
                                                                durations = c(1e-10, 1e10)),
                                    transform = TRUE, re_formula = NA)
    expect_equal(unname(ends[, 1]), draws$up, tolerance = 1e-6, info = m)
    expect_equal(unname(ends[, 2]), draws$low, tolerance = 1e-6, info = m)
    # the table reports each on its own scale (k is the model's own k, not
    # rescaled to another curve's steepness)
    med <- function(pattern) tab$median[grepl(pattern, tab$parameter)]
    expect_equal(med("^low"), stats::median(draws$low), tolerance = 1e-8, info = m)
    expect_equal(med("^up"), stats::median(draws$up), tolerance = 1e-8, info = m)
    expect_equal(med("^k "), stats::median(draws$k), tolerance = 1e-8, info = m)
    # Recovery. In ten seeds of this design low sat 0.002-0.035 below its truth
    # in every matching fit, and its 99.8% interval missed the truth in 4 of 40
    # fits; k sat below its truth in 36 of 40 (its 99.8% interval missed in 1);
    # up's 99.8% interval held in all 40 (phase3_seed_check.R). Centring the
    # priors on the truth removes most of up's error and of k's, but only a
    # quarter to two fifths of low's. So up must lie in its 99.8% interval; low
    # and k within coarse tolerances, with low not pinned at its bound (which
    # the low tolerance alone would let through).
    expect_true(in_wide(draws$up, tr$up), info = paste(m, "up"))
    expect_gt(stats::quantile(draws$low, 0.999), b$low_min + 0.01)
    expect_lt(abs(med("^low") - tr$low), 0.05, label = paste(m, "low"))
    expect_lt(abs(med("^up") - tr$up), 0.05, label = paste(m, "up"))
    expect_lt(abs(log(med("^k ") / tr$k)), log(2), label = paste(m, "k"))
  }
})

test_that("prediction and heat injury run on every model without warnings", {
  skip_unless_brms()
  # A daily cycle inside the assay range (30-36 degC), so no curve is extrapolated.
  trace <- data.frame(time = seq(0, 48, by = 0.5))
  trace$temp <- round(33 + 3 * sin(2 * pi * trace$time / 24), 2)
  for (m in recovery_models) {
    wf <- load_fixture_workflow_model(m)
    expect_no_warning(predict_survival_curves(wf, temps = c(30, 33, 36)))
    expect_no_warning(tls(wf, params = c("z", "ctmax"), t_ref = 60, target_surv = 0.3))
    for (tg in list("relative", 0.3)) {
      expect_no_warning(o <- predict_heat_injury(trace, wf, target_surv = tg, T_c = 29,
                                                 ndraws = 100, seed = 1))
      expect_equal(o$meta$ndraws, 100L)
      for (meth in c("individual", "cumulative")) {
        expect_no_warning(o <- predict_heat_injury(trace, wf, target_surv = tg, T_c = 29,
                                                   ndraws = 100, seed = 1, shape = "varying",
                                                   varying_method = meth))
        expect_equal(o$meta$ndraws, 100L)
        expect_true(all(is.finite(o$summary$surv_median)), info = paste(m, meth))
      }
    }
  }
})

test_that("LOO runs on every model", {
  skip_unless_brms()
  # brms computes the beta-binomial log-likelihood with extraDistr
  skip_if_not_installed("extraDistr")
  for (m in recovery_models) {
    wf <- load_fixture_workflow_model(m)
    # LOO may warn about high Pareto k, which is a diagnostic, not a failure
    lo <- suppressWarnings(brms::loo(get_brmsfit(wf)))
    expect_s3_class(lo, "loo")
    expect_true(is.finite(lo$estimates["elpd_loo", "Estimate"]), info = m)
  }
})

test_that("a wrong curve barely moves z or CTmax when the design pins the curve", {
  skip_unless_brms()
  # Paired: the same data fitted with the matching curve and with a wrong one,
  # in both directions. Every model shares the midpoint line, so in a design
  # with durations on each curve's shoulders a wrong curve shifts the fitted
  # midpoints by nearly the same amount at every temperature. Ten seeds
  # (phase3_seed_check.R): the paired shifts were at most 0.24 in z and
  # 0.40 degC in CTmax. At 63 rows the finite-sample and prior effects are as
  # large as the curve's own (exact large-sample shifts for these pairs are
  # within 0.14 in z and 0.20 degC in CTmax; kl_design_check.R). In a sparse
  # design (decade steps) the data fix z for no curve at all, so the check
  # needs a design like this one.
  for (pr in list(c("weibull1", "log-logistic"), c("log-logistic", "weibull1"))) {
    match <- load_fixture_workflow_model(pr[1])
    wrong <- load_fixture_workflow_model(pr[1], fit_model = pr[2])
    expect_equal(wrong$meta$model, pr[2])
    tr  <- attr(wrong, "truth")
    sm  <- tls(match, params = c("z", "ctmax"), t_ref = 60)
    sw  <- tls(wrong, params = c("z", "ctmax"), t_ref = 60)
    med <- function(s, q) s$summary$median[s$summary$quantity == q]
    lab <- paste(pr[1], "data,", pr[2], "fit")
    expect_lt(abs(med(sw, "z") - med(sm, "z")), 0.5, label = paste(lab, "z"))
    expect_lt(abs(med(sw, "CTmax") - med(sm, "CTmax")), 0.6, label = paste(lab, "CTmax"))
    expect_true(in_wide(sw$draws$value[sw$draws$quantity == "z"], tr$z), info = lab)
  }
})
