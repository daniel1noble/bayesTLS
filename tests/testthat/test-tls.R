# Accuracy tests for tls(). The brms-free block checks argument validation; the
# gated block validates numerical correctness against the exact -1/beta1
# identity and the simulation truth (z, CTmax, T_crit), plus the per-group
# machinery, threshold edge cases and ndraws clamping.

test_that("tls() validates arguments before touching a fit", {
  expect_error(tls(list(foo = 1)), "bayes_tls|brmsfit")
  expect_error(tls(list(foo = 1), params = "tcrit", lethal = FALSE), "lethal")
})

test_that("tls() relative z equals the exact -1/b_mid_temp_c identity and recovers truth", {
  skip_unless_brms()
  wf <- load_fixture_workflow()
  tr <- truth_summary()
  d  <- posterior::as_draws_df(get_brmsfit(wf))

  tl    <- tls(wf, params = "z", mode = "relative")
  z_tls <- tl$summary$median[tl$summary$quantity == "z"]

  # Exact closed form: relative z = -1 / b_mid_temp_c per draw (tls uses all draws).
  expect_equal(z_tls, stats::median(-1 / d$b_mid_temp_c), tolerance = 1e-6)
  # Recovers the known simulation truth (z = -1/m_beta1 = 5.556).
  expect_equal(z_tls, tr$z, tolerance = 0.5)
})

test_that("tls() recovers z and CTmax_1hr from simulated truth", {
  skip_unless_brms()
  wf  <- load_fixture_workflow()
  tr  <- truth_summary()
  out <- tls(wf, params = c("z", "ctmax"), t_ref = 60, ndraws = 500)$summary
  g   <- function(q, col) out[[col]][out$quantity == q]

  # Posterior medians within tolerance of truth ...
  expect_lt(abs(g("z", "median")     - tr$z),         1.5)
  expect_lt(abs(g("CTmax", "median") - tr$CTmax_1hr), 1.0)
  # ... and truth inside the 95% credible interval.
  expect_lte(g("z", "lower"),     tr$z)
  expect_gte(g("z", "upper"),     tr$z)
  expect_lte(g("CTmax", "lower"), tr$CTmax_1hr)
  expect_gte(g("CTmax", "upper"), tr$CTmax_1hr)
})

test_that("tls(lethal = TRUE) returns a sensible T_crit below CTmax", {
  skip_unless_brms()
  wf  <- load_fixture_workflow()
  tr  <- truth_summary()
  out <- suppressMessages(
    tls(wf, params = "all", t_ref = 60, ndraws = 500, lethal = TRUE))$summary
  g   <- function(q, col) out[[col]][out$quantity == q]

  # With default TC_rate_range = c(0.1, 1) the T_crit median sits near
  # CTmax_1hr - 2.5 * z (geometric mean of the log10(rate) range); the true
  # T_crit median must fall inside the posterior 95% CrI, and below CTmax.
  true_T_crit_median <- tr$CTmax_1hr - 2.5 * tr$z
  expect_lt(g("Tcrit", "median"), g("CTmax", "median"))
  expect_lte(g("Tcrit", "lower"), true_T_crit_median)
  expect_gte(g("Tcrit", "upper"), true_T_crit_median)
  # CrI width > 0.5 z (pooling over log10(r*) contributes ~1 z of spread).
  expect_gt(g("Tcrit", "upper") - g("Tcrit", "lower"), 0.5 * tr$z)
})

test_that("an absolute target_surv outside the fitted asymptotes warns clearly", {
  skip_unless_brms()
  wf <- load_fixture_workflow_grouped()   # survival spans ~[0.04, 0.97]
  # 0.1 lies below the lower asymptote for some draws -> LT undefined there:
  # a clear, reasoned warning instead of a bare "NaNs produced".
  expect_warning(
    muffle_uncrossed(tls(wf, by = "grp", target_surv = 0.1, ndraws = NULL,
                         temp_mean = wf$meta$temp_mean)),
    "outside the fitted curve's asymptotes")
  # an in-range absolute threshold does not raise that warning
  asymp_warn <- FALSE
  withCallingHandlers(
    tls(wf, by = "grp", target_surv = 0.5, ndraws = NULL,
        temp_mean = wf$meta$temp_mean),
    warning = function(w) {
      if (grepl("asymptotes", conditionMessage(w))) asymp_warn <<- TRUE
      invokeRestart("muffleWarning")
    })
  expect_false(asymp_warn)
})

test_that("derive_tdt_curve() clamps ndraws to the posterior size", {
  skip_unless_brms()
  wf <- load_fixture_workflow()                 # ~1500 posterior draws
  # Regression: relative-mode derive_tdt_curve() once passed ndraws straight to
  # brms::posterior_linpred(), which errors when ndraws exceeds the posterior
  # size. Requesting more draws than exist must clamp, not crash.
  expect_no_error(derive_tdt_curve(wf, temp_grid = c(30, 34, 38), ndraws = 5000))
})

test_that("tls() params/lethal switches and summary shape", {
  skip_unless_brms()
  wf <- load_fixture_workflow()

  z_only <- tls(wf, params = "z", mode = "relative")
  expect_s3_class(z_only, "tls")
  expect_equal(unique(z_only$summary$quantity), "z")
  expect_named(z_only$summary, c("quantity", "median", "lower", "upper"))

  all3 <- tls(wf, params = "all", lethal = TRUE, mode = "relative")
  expect_setequal(unique(all3$summary$quantity), c("z", "CTmax", "Tcrit"))
  # draws slot carries one row per draw per quantity
  expect_true(all(c("quantity", ".draw", "value") %in% names(all3$draws)))
})

test_that("tls() warns and NAs the summary when a group's z/CTmax draws are all non-finite", {
  skip_unless_brms()
  wf  <- load_fixture_workflow()
  fit <- get_brmsfit(wf)
  # A single-temperature grid makes the LS slope 0/0 (NaN) -> z = -1/NaN non-finite
  # for every draw, so the group's summary must be NA and a naming warning fired.
  nd  <- data.frame(grp = factor("solo"),
                    temp_c = mean(fit$data$temp_c))

  expect_warning(
    out <- tls(wf, by = "grp", newdata = nd, params = "z", mode = "relative"),
    "non-finite.*solo|solo.*non-finite"
  )
  expect_true(all(is.na(out$summary$median[out$summary$quantity == "z"])))
  # Raw draws are still returned (non-finite, not silently dropped).
  expect_true(all(!is.finite(out$draws$value[out$draws$quantity == "z"])))
})

test_that("tls() per-group machinery: a moderator the model ignores gives identical z", {
  skip_unless_brms()
  wf  <- load_fixture_workflow()
  fit <- get_brmsfit(wf)
  tc  <- seq(min(fit$data$temp_c), max(fit$data$temp_c), length.out = 5)
  # `grp` is NOT in the model, so both groups must yield the same dose-response.
  nd  <- expand.grid(grp = factor(c("a", "b")), temp_c = tc)

  byg <- tls(wf, by = "grp", newdata = nd, params = "z", mode = "relative")$summary
  za  <- byg$median[byg$grp == "a" & byg$quantity == "z"]
  zb  <- byg$median[byg$grp == "b" & byg$quantity == "z"]
  pooled <- tls(wf, params = "z", mode = "relative")$summary$median

  expect_equal(za, zb, tolerance = 1e-8)       # model ignores grp -> identical
  expect_equal(za, pooled, tolerance = 0.02)   # and equals the pooled estimate
})

# --- by = NULL uses the fit's own groups (fixed 2026-11-26) -------------------
# Previously tls() with by = NULL on a grouped fit filled the moderator from
# data row 1 and returned that one group under an ungrouped label.

test_that("tls_resolve_by(): explicit by wins, newdata is taken as given, workflows use their groups", {
  grouped <- structure(list(fit = NULL, meta = list(group_vars = "grp", grouped = TRUE)),
                       class = "bayes_tls")
  single  <- structure(list(fit = NULL, meta = list(group_vars = character(0),
                                                    grouped = FALSE)),
                       class = "bayes_tls")
  expect_equal(tls_resolve_by(grouped, NULL, NULL), "grp")
  expect_equal(tls_resolve_by(grouped, NULL, "other"), "other")
  expect_null(tls_resolve_by(grouped, NULL, NULL, newdata = data.frame(temp_c = 0)))
  expect_null(tls_resolve_by(single, NULL, NULL))
  # A grouped-looking workflow that records no moderator must name it.
  unnamed <- structure(list(fit = NULL, meta = list(grouped = TRUE)), class = "bayes_tls")
  expect_error(tls_resolve_by(unnamed, NULL, NULL), "Pass `by")
})

test_that("tls() without by reports every group of a grouped workflow (gated)", {
  skip_unless_brms()
  for (wf in list(load_fixture_workflow_2group_beta_direct(),
                  load_fixture_workflow_2group_binom_mid())) {
    auto <- tls(wf, params = c("z", "ctmax"), t_ref = 60, ndraws = 300, seed = 3)$summary
    expl <- tls(wf, by = "grp", params = c("z", "ctmax"), t_ref = 60, ndraws = 300,
                seed = 3)$summary
    expect_true("grp" %in% names(auto))
    expect_equal(nrow(auto), 4)
    expect_equal(auto, expl)
  }
})

test_that("tls() on a bare grouped brmsfit asks for by (gated)", {
  skip_unless_brms()
  wf  <- load_fixture_workflow_2group_beta_direct()
  fit <- get_brmsfit(wf)
  expect_error(tls(fit, params = "z", temp_mean = wf$meta$temp_mean), "Pass `by")
  s <- tls(fit, by = "grp", params = "z", temp_mean = wf$meta$temp_mean, t_ref = 60)$summary
  expect_equal(nrow(s), 2)
})

# --- Exact absolute crossings and shared draw subsamples (fixed 2026-11-26) ----

test_that("tls_check_temp_grid() warns only for a grid far outside the data", {
  obs <- c(-3, 0, 3)                                   # centred temperatures
  expect_silent(tls_check_temp_grid(NULL, obs, 33))
  expect_silent(tls_check_temp_grid(c(-3, 0, 3), obs, 33))
  expect_silent(tls_check_temp_grid(c(-8, 8), obs, 33))  # within one range width
  expect_warning(tls_check_temp_grid(c(30, 33, 36), obs, 33), "temp_grid - 33")
  expect_warning(tls_check_temp_grid(c(-12, 0), obs, NULL), "far outside")
})

test_that("tls() warns when temp_grid looks like degrees C (gated)", {
  skip_unless_brms()
  wf <- load_fixture_workflow_beta()
  expect_warning(tls(wf, params = "z", t_ref = 60, temp_grid = c(30, 33, 36)),
                 "centring temperature")
  expect_no_warning(tls(wf, params = "z", t_ref = 60,
                        temp_grid = c(30, 33, 36) - wf$meta$temp_mean))
})

test_that("tls() uses one draw subsample for all four sub-parameters (gated)", {
  skip_unless_brms()
  wf  <- load_fixture_workflow_beta()
  fit <- get_brmsfit(wf)
  sub  <- muffle_uncrossed(tls(wf, params = c("z", "ctmax"), t_ref = 60,
                               target_surv = 0.3, ndraws = 60, seed = 4))$draws
  full <- muffle_uncrossed(tls(wf, params = c("z", "ctmax"), t_ref = 60,
                               target_surv = 0.3))$draws
  set.seed(4)
  ids <- tls_draw_ids(fit, 60)
  for (q in c("z", "CTmax")) {
    a <- sub[sub$quantity == q, ]
    b <- full[full$quantity == q, ]
    expect_equal(a$value[order(a$.draw)], b$value[match(ids, b$.draw)])
  }
})

test_that("derive_tdt_curve() reads absolute crossings exactly, beyond any grid (gated)", {
  skip_unless_brms()
  wf  <- load_fixture_workflow_beta()
  fit <- get_brmsfit(wf)
  temps <- c(15, 33)                    # 15 C crosses far beyond the data's durations
  nd  <- new_tdt_grid(wf, temps = temps, durations = 1)
  for (p in c(0.3, 0.5)) {              # 0.5 is also the default: test another p too
    # At p = 0.3 some 15 C cells lie outside the asymptotes (the engine warns)
    cur <- suppressWarnings(derive_tdt_curve(wf, temp_grid = temps, target_surv = p,
                                             time_multiplier = 1,
                                             output_time_unit = "hours"))
    eng <- suppressWarnings(tls_eval_subpars(fit, nd, wf$meta$bounds,
                                             mode = "absolute", p = p))$logLT
    for (j in seq_along(temps)) {
      d <- cur$draws[cur$draws$temp == temps[j], ]
      expect_equal(nrow(d), sum(is.finite(eng[, j])))        # only undefined draws dropped
      expect_equal(d$duration_model, 10^eng[d$.draw, j], tolerance = 1e-10)
      if (p == 0.5) expect_equal(nrow(d), brms::ndraws(fit)) # none at p = 0.5
    }
  }
  cur <- derive_tdt_curve(wf, temp_grid = temps, target_surv = 0.5,
                          time_multiplier = 1, output_time_unit = "hours")
  expect_gt(cur$summary$duration_median[cur$summary$temp == 15],
            5 * max(wf$data$duration))                       # beyond the old grid
  # With ndraws, the four sub-parameters share one subsample.
  set.seed(9)
  ids <- tls_draw_ids(fit, 40)
  set.seed(9)
  cur40 <- derive_tdt_curve(wf, temp_grid = temps, target_surv = 0.5, ndraws = 40,
                            time_multiplier = 1, output_time_unit = "hours")
  eng40 <- tls_eval_subpars(fit, nd, wf$meta$bounds, draw_ids = ids,
                            mode = "absolute", p = 0.5)$logLT
  d40 <- cur40$draws[cur40$draws$temp == 33, ]
  expect_equal(log10(d40$duration_model), eng40[d40$.draw, 2], tolerance = 1e-10)
  # The retired search grid is ignored, with a message, in either mode.
  expect_message(derive_tdt_curve(wf, temp_grid = 33, target_surv = 0.5,
                                  duration_grid = c(1, 10)), "no longer used")
  expect_message(derive_tdt_curve(wf, temp_grid = 33, duration_grid = c(1, 10)),
                 "no longer used")
})

test_that("derive_tdt_curve() absolute is wired per group (gated)", {
  skip_unless_brms()
  wf  <- load_fixture_workflow_2group_beta_direct()
  fit <- get_brmsfit(wf)
  cur <- derive_tdt_curve(wf, temp_grid = c(30, 36), target_surv = 0.3, by = "grp",
                          time_multiplier = 1, output_time_unit = "hours")
  nd  <- new_tdt_grid(wf, temps = c(30, 36), durations = 1, by = "grp")
  eng <- tls_eval_subpars(fit, nd, wf$meta$bounds, mode = "absolute", p = 0.3)$logLT
  expect_setequal(unique(cur$draws$grp), unique(as.character(nd$grp)))
  for (i in seq_len(nrow(nd))) {
    d <- cur$draws[cur$draws$grp == nd$grp[i] & cur$draws$temp == nd$temp[i], ]
    expect_equal(log10(d$duration_model), eng[d$.draw, i], tolerance = 1e-10)
  }
  # the two groups really differ at 36 C (so a collapsed grouping would show)
  s36 <- cur$summary[cur$summary$temp == 36, ]
  expect_gt(abs(diff(log10(s36$duration_median))), 0.1)
})

# --- CTmax crossings outside the temperature grid (fixed 2026-11-26) ----------

test_that("tls() reads a direct fit's relative CTmax in closed form, with no draw dropped (gated)", {
  skip_unless_brms()
  wd <- load_fixture_workflow_direct()
  fd <- get_brmsfit(wd)
  mid <- brms::posterior_linpred(fd, newdata = tls_build_grid(fd$data, temp_grid = c(0, 1)),
                                 nlpar = "mid", re_formula = NA)
  tm  <- wd$meta$temp_mean
  for (tr in c(6, 60)) {                                 # 6 min: CTmax beyond the data
    v  <- tls(wd, params = "ctmax", t_ref = tr)$draws$value
    cf <- tm + (log10(tr / 60) - mid[, 1]) / (mid[, 2] - mid[, 1])
    expect_true(all(is.finite(v)))
    expect_equal(v, cf, tolerance = 1e-9)
    if (tr == 6) expect_gt(max(cf), max(wd$data$temp))   # some draws lie off the grid
  }
  tt <- derive_temperature_for_duration(wd, exposure_duration = 1)$draws$temp
  expect_equal(tt, tls(wd, params = "ctmax", t_ref = 60)$draws$value, tolerance = 1e-9)
})

test_that("tls() inverts a bent curve finely and searches outward, matching a brute-force search (gated)", {
  skip_unless_brms()
  wb  <- load_fixture_workflow_beta()
  fit <- get_brmsfit(wb); tm <- wb$meta$temp_mean
  target <- log10(6 / 60)                                # t_ref = 6 min, model in hours
  expect_warning(res <- tls(wb, params = "ctmax", target_surv = 0.3, t_ref = 6),
                 "no single crossing of their target exposure")
  new <- res$draws$value
  tc   <- fit$data$temp_c
  w    <- diff(range(tc))                                # the outward search's reach
  fine <- seq(min(tc) - w, max(tc) + w, by = 0.005)
  sp <- suppressWarnings(tls_eval_subpars(fit, tls_build_grid(fit$data, temp_grid = fine),
                                          wb$meta$bounds, mode = "absolute", p = 0.3))
  Mf  <- sp$logLT
  ins <- which(fine >= min(tc) - 1e-9 & fine <= max(tc) + 1e-9)
  ref <- vapply(seq_len(nrow(Mf)), function(i) {   # single falling crossing in the data's
    y  <- Mf[i, ins] - target; n <- length(y)       # range, else the first one outward
    j  <- which(y[-n] > 0 & y[-1] <= 0)
    if (length(j) == 1) { t1 <- fine[ins][j]; t2 <- fine[ins][j + 1]
      return(tm + t1 + y[j] * (t2 - t1) / (y[j] - y[j + 1])) }
    if (length(j) > 1 || anyNA(y)) return(NA_real_)
    cols <- if (all(y > 0)) which(fine >= max(tc) - 1e-9) else
      if (all(y <= 0)) rev(which(fine <= min(tc) + 1e-9)) else return(NA_real_)
    yy <- Mf[i, cols] - target; k <- which(is.na(yy) | (if (all(y > 0)) yy <= 0 else yy > 0))[1]
    if (is.na(k) || is.na(yy[k])) return(NA_real_)
    tt <- fine[cols]; tm + tt[k - 1] + yy[k - 1] * (tt[k] - tt[k - 1]) / (yy[k - 1] - yy[k])
  }, numeric(1))
  expect_gt(sum(is.finite(ref)), 1400)                   # most draws cross, many off the grid
  expect_lt(mean(xor(is.finite(new), is.finite(ref))), 0.005)
  ok <- is.finite(new) & is.finite(ref)
  # 3e-4 on this (git-tracked) fixture; a coarse grid or a coarse outward step
  # gives 1e-3 to 3e-3
  expect_lt(max(abs(new[ok] - ref[ok])), 1e-3)
  # crossings inside the data are read on the refined default grid
  tf  <- tls_refine_grid(seq(min(tc), max(tc), length.out = 11))
  Mi  <- suppressWarnings(tls_eval_subpars(fit, tls_build_grid(fit$data, temp_grid = tf),
                                           wb$meta$bounds, mode = "absolute", p = 0.3))$logLT
  inv <- tls_invert_logLT(Mi, target, tf + tm)
  ins <- is.finite(inv)
  expect_gt(sum(ins), 0)
  expect_equal(new[ins], inv[ins], tolerance = 1e-10)
})

test_that("tls() and derive_temperature_for_duration() agree per draw on a bent curve (gated)", {
  skip_unless_brms()
  wb <- load_fixture_workflow_beta()
  a <- suppressWarnings(tls(wb, params = "ctmax", target_surv = 0.3, t_ref = 6))$draws$value
  expect_warning(b <- derive_temperature_for_duration(wb, exposure_duration = 6 / 60,
                                                      target_surv = 0.3),
                 "no single crossing of the exposure")
  bb <- rep(NA_real_, length(a)); bb[b$draws$.draw] <- b$draws$temp
  expect_equal(bb, a, tolerance = 1e-8)
})

test_that("tls() with ndraws keeps the outward search on the right draws (gated)", {
  skip_unless_brms()
  wb  <- load_fixture_workflow_beta()
  fit <- get_brmsfit(wb)
  sub  <- suppressWarnings(tls(wb, params = "ctmax", target_surv = 0.3, t_ref = 6,
                               ndraws = 80, seed = 11))$draws
  full <- suppressWarnings(tls(wb, params = "ctmax", target_surv = 0.3, t_ref = 6))$draws
  set.seed(11)
  ids <- tls_draw_ids(fit, 80)
  expect_equal(sub$value, full$value[ids])
})

test_that("tls() and the readers reject an invalid ndraws (gated)", {
  skip_unless_brms()
  wf <- load_fixture_workflow_beta()
  for (bad in list(0, -5, NA, 2.5)) {
    expect_error(tls(wf, params = "z", t_ref = 60, ndraws = bad), "whole number")
    expect_error(derive_tdt_curve(wf, temp_grid = 33, ndraws = bad), "whole number")
    expect_error(derive_tdt_curve(wf, temp_grid = 33, target_surv = 0.5, ndraws = bad),
                 "whole number")
  }
  big <- tls(wf, params = "z", t_ref = 60, ndraws = 1e6)$draws    # more than the posterior
  expect_equal(nrow(big), brms::ndraws(get_brmsfit(wf)))
})

test_that("readers work on a hand-built workflow without meta$bounds (gated)", {
  skip_unless_brms()
  wf <- load_fixture_workflow_beta()
  wrap <- wf
  wrap$meta$bounds <- NULL                               # e.g. a hand-built bayes_tls
  expect_equal(derive_tdt_curve(wrap, temp_grid = c(30, 36), target_surv = 0.5)$summary,
               derive_tdt_curve(wf,   temp_grid = c(30, 36), target_surv = 0.5)$summary)
  expect_equal(derive_z(wrap, target_surv = 0.5, seed = 1)$summary,
               derive_z(wf,   target_surv = 0.5, seed = 1)$summary)
  expect_equal(tdt_parameter_table(wrap), tdt_parameter_table(wf))
})

# --- Absolute targets near an asymptote; search wiring (second review) --------

# Reference CTmax per draw, independent of the LT: survival at t_ref minus p on
# a fine temperature grid (defined everywhere, also where the LT at p is not),
# taking the single falling crossing inside the data's temperatures, else the
# first crossing outward (up to one data range beyond each end).
ref_surv_ctmax <- function(wf, p, t_ref, step = 0.01, by = NULL) {
  fit <- get_brmsfit(wf); tm <- wf$meta$temp_mean
  target <- log10(t_ref / tdt_resolve_time_multiplier(NULL, wf$meta, "min"))
  tc <- fit$data$temp_c; lo <- min(tc); hi <- max(tc); w <- hi - lo
  fine <- seq(lo - w, hi + w, by = step)
  nd <- tls_build_grid(fit$data, by = by, temp_grid = fine)
  sp <- suppressWarnings(tls_eval_subpars(fit, nd, wf$meta$bounds, mode = "absolute", p = p))
  rr <- tdt_as_model_spec(tdt_model(fit))$curve(sp$low, sp$up, sp$k * (target - sp$mid)) - p
  ins <- which(fine >= lo - 1e-9 & fine <= hi + 1e-9)
  out <- lapply(split(seq_len(nrow(nd)), nd$.grp), function(gi) apply(rr[, gi, drop = FALSE], 1, function(y) {
    yi <- y[ins]; j <- which(yi[-length(yi)] > 0 & yi[-1] <= 0)
    if (length(j) == 1) { t1 <- fine[ins][j]; t2 <- fine[ins][j + 1]
      return(tm + t1 + yi[j] * (t2 - t1) / (yi[j] - yi[j + 1])) }
    if (length(j) > 1) return(NA_real_)
    up_ <- all(yi > 0)
    cols <- if (up_) which(fine >= hi - 1e-9) else if (all(yi <= 0)) rev(which(fine <= lo + 1e-9)) else return(NA_real_)
    yy <- y[cols]; k <- which(if (up_) yy <= 0 else yy > 0)[1]
    if (is.na(k)) return(NA_real_)
    t1 <- fine[cols][k - 1]; t2 <- fine[cols][k]
    tm + t1 + yy[k - 1] * (t2 - t1) / (yy[k - 1] - yy[k])
  }))
  if (is.null(by)) out[[1]] else out
}

test_that("tls() finds CTmax next to an asymptote, matching a survival-curve reference (gated)", {
  skip_unless_brms()
  ws <- load_fixture_workflow_shape()          # temperature on all four parameters
  for (p in c(0.8, 0.85)) {
    got <- suppressWarnings(tls(ws, params = "ctmax", target_surv = p, t_ref = 60))$draws$value
    ref <- ref_surv_ctmax(ws, p, 60)
    expect_identical(is.finite(got), is.finite(ref))
    expect_equal(sum(is.finite(got)), length(got))      # every draw has a CTmax (2 at HEAD at 0.8)
    expect_lt(max(abs(got - ref)), 0.1)                  # refined grid next to the asymptote
    # 0.006-0.008 on fixtures refitted with other seeds (about 2% of refits
    # would exceed 0.01, from draws crossing next to the asymptote)
    expect_lt(max(abs(stats::quantile(got, c(.025, .5, .975)) -
                      stats::quantile(ref, c(.025, .5, .975)))), 0.02)
    dd <- suppressWarnings(derive_temperature_for_duration(ws, 1, target_surv = p))$draws
    expect_equal(dd$temp, got[dd$.draw], tolerance = 1e-8)
  }
  wb  <- load_fixture_workflow_beta()           # upper asymptote below 0.9 at the hot end
  got <- muffle_uncrossed(suppressWarnings(tls(wb, params = "ctmax", target_surv = 0.9, t_ref = 60)))
  v   <- got$draws$value
  ref <- ref_surv_ctmax(wb, 0.9, 60)
  expect_gt(sum(is.finite(ref)), 500)                   # most draws cross (636 of 1500 kept before this fix)
  expect_lt(mean(xor(is.finite(v), is.finite(ref))), 0.01)
  ok <- is.finite(v) & is.finite(ref)
  expect_lt(max(abs(v[ok] - ref[ok])), 0.1)             # crossings next to the asymptote
})

test_that("T_crit's 1-hour anchor gets the same outward search, and the warning counts per quantity (gated)", {
  skip_unless_brms()
  wb <- load_fixture_workflow_beta()
  # TC_rate_range = c(100, 100) makes the rate term zero, so T_crit = its anchor
  a <- muffle_uncrossed(tls(wb, params = "tcrit", lethal = TRUE, target_surv = 0.3,
                            t_ref = 6, TC_rate_range = c(100, 100)))$draws$value
  d <- muffle_uncrossed(derive_temperature_for_duration(wb, 1, target_surv = 0.3))$draws
  full <- rep(NA_real_, length(a)); full[d$.draw] <- d$temp
  expect_equal(a, full, tolerance = 1e-8)
  ref <- ref_surv_ctmax(wb, 0.3, 60)                    # the anchor is the 1-hour CTmax
  expect_lt(mean(xor(is.finite(a), is.finite(ref))), 0.01)
  w <- character(0)
  r <- withCallingHandlers(
    tls(wb, params = c("ctmax", "tcrit"), lethal = TRUE, target_surv = 0.3, t_ref = 20, seed = 1),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  nc <- sum(!is.finite(r$draws$value[r$draws$quantity == "CTmax"]))
  nt <- sum(!is.finite(r$draws$value[r$draws$quantity == "Tcrit"]))
  uw <- grep("no single crossing", w, value = TRUE)
  if (nc + nt == 0) {
    expect_length(uw, 0L)
  } else {                                              # one warning, counted per quantity
    expect_length(uw, 1L)
    if (nc > 0) expect_match(uw, sprintf("CTmax %d of 1500 draws", nc), fixed = TRUE)
    else expect_false(grepl("CTmax [0-9]+ of", uw))
    if (nt > 0) expect_match(uw, sprintf("T_crit %d of 1500 draws", nt), fixed = TRUE)
    else expect_false(grepl("T_crit [0-9]+ of", uw))
    expect_match(uw, "searched 24 to 42", fixed = TRUE)
  }
})

test_that("the outward search reaches the data's range from a narrow grid, and derive_* is wired alike (gated)", {
  skip_unless_brms()
  wb  <- load_fixture_workflow_beta()
  fit <- get_brmsfit(wb)
  narrow <- muffle_uncrossed(tls(wb, params = "ctmax", target_surv = 0.3, t_ref = 6,
                                 temp_grid = c(-1, 0, 1)))$draws$value
  wide   <- muffle_uncrossed(tls(wb, params = "ctmax", target_surv = 0.3,
                                 t_ref = 6))$draws$value
  tm <- wb$meta$temp_mean                               # data 30-36 C: reach -1..1 +/- 6
  reach <- is.finite(wide) & wide > tm - 7 + 0.1 & wide < tm + 7 - 0.1
  expect_gt(sum(reach & wide > tm + 3), 100)            # beyond the grid's own width of 2 C
  expect_gt(mean(is.finite(narrow[reach])), 0.99)
  expect_true(all(abs(narrow - tm) <= 7 + 1e-8, na.rm = TRUE))
  dn <- muffle_uncrossed(derive_temperature_for_duration(wb, 0.1, target_surv = 0.3,
                                                         temp_grid = c(32, 33, 34)))$draws
  expect_equal(dn$temp, narrow[dn$.draw], tolerance = 1e-8)
  expect_warning(derive_temperature_for_duration(wb, 0.1, target_surv = 0.3,
                                                 temp_grid = c(32, 33, 34)),
                 "within 26 to 40 degrees C", fixed = TRUE)   # -1..1 +/- 6 (data range)
  # a reversed grid gives the same result; ndraws keeps the right draws
  rev_g <- muffle_uncrossed(derive_temperature_for_duration(
    wb, 0.1, target_surv = 0.3, temp_grid = rev(seq(30, 36, length.out = 11))))$draws
  fwd_g <- muffle_uncrossed(derive_temperature_for_duration(
    wb, 0.1, target_surv = 0.3, temp_grid = seq(30, 36, length.out = 11)))$draws
  expect_equal(rev_g, fwd_g)
  sub  <- muffle_uncrossed(derive_temperature_for_duration(wb, 0.1, target_surv = 0.3,
                                                           ndraws = 80, seed = 11))$draws
  full <- muffle_uncrossed(derive_temperature_for_duration(wb, 0.1, target_surv = 0.3))$draws
  set.seed(11); ids <- tls_draw_ids(fit, 80)
  fv <- rep(NA_real_, brms::ndraws(fit)); fv[full$.draw] <- full$temp
  sv <- rep(NA_real_, 80); sv[sub$.draw] <- sub$temp
  expect_equal(sv, fv[ids])
})

test_that("tls() drops non-finite temp_grid values with a warning (gated)", {
  skip_unless_brms()
  wf <- load_fixture_workflow()
  expect_warning(a <- tls(wf, params = "ctmax", t_ref = 60, temp_grid = c(-3, 0, 3, Inf)),
                 "non-finite")
  expect_equal(a$draws, tls(wf, params = "ctmax", t_ref = 60, temp_grid = c(-3, 0, 3))$draws)
})

test_that("more readers work on a workflow without meta$bounds (gated)", {
  skip_unless_brms()
  wf <- load_fixture_workflow_beta()
  wrap <- wf
  wrap$meta$bounds <- NULL
  expect_equal(muffle_uncrossed(derive_temperature_for_duration(wrap, 1, target_surv = 0.5))$summary,
               muffle_uncrossed(derive_temperature_for_duration(wf,   1, target_surv = 0.5))$summary)
  expect_equal(extract_4pl_pars(wrap), extract_4pl_pars(wf))
  tr <- data.frame(time = seq(0, 2, by = 0.25), temp = 33)
  expect_equal(predict_heat_injury(tr, wrap, ndraws = 50, seed = 1)$summary,
               predict_heat_injury(tr, wf,   ndraws = 50, seed = 1)$summary)
})

test_that("the batched outward search keeps each group's own extension (gated)", {
  skip_unless_brms()
  wg <- load_fixture_workflow_shape_grouped()   # both groups cross above the data at p = 0.3
  # (p, t_ref): at 0.3 / 120 min group A needs the outward search for all 1000
  # draws and group B for a subset of 814, so the shared batch must map rows
  # to each group's own draws.
  for (cs in list(c(0.3, 60), c(0.3, 120), c(0.8, 60))) {
    p <- cs[1]; tr <- cs[2]
    r   <- muffle_uncrossed(suppressWarnings(tls(wg, params = "ctmax", target_surv = p,
                                                 t_ref = tr)))$draws
    ref <- ref_surv_ctmax(wg, p, tr, by = "grp")
    for (g in names(ref)) {
      got <- r$value[as.character(r$grp) == g]
      lab <- paste("group", g, "p", p, "t_ref", tr)
      expect_lt(mean(xor(is.finite(got), is.finite(ref[[g]]))), 0.01, label = lab)
      ok  <- is.finite(got) & is.finite(ref[[g]])
      expect_gt(sum(ok), 500, label = lab)
      expect_lt(max(abs(got[ok] - ref[[g]][ok])), 0.1, label = lab)    # refined grid
      expect_lt(abs(stats::median(got[ok]) - stats::median(ref[[g]][ok])), 0.01, label = lab)
    }
  }
})

test_that("derive_temperature_for_duration() inverts grouped bent curves like tls() (gated)", {
  skip_unless_brms()
  wg  <- load_fixture_workflow_shape_grouped()
  d   <- muffle_uncrossed(suppressWarnings(
    derive_temperature_for_duration(wg, 1, target_surv = 0.8, by = "grp")))$draws
  r   <- muffle_uncrossed(suppressWarnings(tls(wg, params = "ctmax", target_surv = 0.8,
                                               t_ref = 60)))$draws
  ref <- ref_surv_ctmax(wg, 0.8, 60, by = "grp")
  for (g in names(ref)) {
    got <- r$value[as.character(r$grp) == g]
    dg  <- d[as.character(d$grp) == g, ]
    expect_equal(dg$temp, got[dg$.draw], tolerance = 1e-8)
    expect_equal(sum(is.finite(got)), nrow(dg))
    ok <- is.finite(got) & is.finite(ref[[g]])
    expect_lt(mean(xor(is.finite(got), is.finite(ref[[g]]))), 0.01, label = g)
    expect_lt(max(abs(got[ok] - ref[[g]][ok])), 0.1, label = g)
  }
  # no crossing in any group: an NA row per group and the warning, not an error
  w <- character(0)
  none <- withCallingHandlers(
    derive_temperature_for_duration(wg, 1e8, target_surv = 0.5, by = "grp"),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("no single crossing", w, fixed = TRUE)))
  expect_equal(nrow(none$draws), 0L)
  expect_setequal(as.character(none$summary$grp), c("A", "B"))
  expect_true(all(is.na(none$summary$temp_median)))
})

test_that("derive_temperature_for_duration() is not fooled by a grid at temp_mean and temp_mean + 1 (gated)", {
  skip_unless_brms()
  ws  <- load_fixture_workflow_shape()          # bent: temperature on all four parameters
  tm  <- ws$meta$temp_mean
  ref <- ref_surv_ctmax(ws, 0.8, 600)           # 10 h: crossings inside the data
  for (g in list(tm, tm + 1, c(tm, tm + 1))) {
    d <- muffle_uncrossed(suppressWarnings(
      derive_temperature_for_duration(ws, 10, target_surv = 0.8, temp_grid = g)))$draws
    got <- rep(NA_real_, length(ref)); got[d$.draw] <- d$temp
    ok  <- is.finite(got) & is.finite(ref)
    lab <- paste("temp_grid", paste(g, collapse = ", "))
    expect_gt(sum(ok), 0.95 * length(ref), label = lab)
    expect_lt(max(abs(got[ok] - ref[ok])), 0.05, label = lab)   # the old straight line: 0.16
  }
})

test_that("a group with no CTmax at all gets one warning, not two (gated)", {
  skip_unless_brms()
  ws <- load_fixture_workflow_shape()          # CTmax at 6 min lies beyond the search
  w <- character(0)
  r <- withCallingHandlers(tls(ws, params = "ctmax", target_surv = 0.5, t_ref = 6),
         warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  expect_equal(sum(is.finite(r$draws$value)), 0)
  expect_length(w, 1L)
  expect_match(w, "all z/CTmax draws were non-finite", fixed = TRUE)
  expect_match(w, "within the temperatures searched", fixed = TRUE)
})
