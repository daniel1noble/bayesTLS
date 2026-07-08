# PILLAR 3 -- two interacting CONTINUOUS moderators (x1, x2), centred at 0.
#
# `by=` is factors-only, so the moderators are specified directly:
#   fit_4pl(..., ctmax = ~ x1 * x2, z = ~ x1 * x2).
# The data are simulated from KNOWN direct-parameterisation surfaces
#   CTmaxdev(x1,x2) = a0 + a1 x1 + a2 x2 + a3 x1:x2
#   logz(x1,x2)     = z0 + z1 x1 + z2 x2 + z3 x1:x2
# (see simulate_tdt_cont2()). Recovery is asserted three ways:
#   (a) COEFFICIENT level    -- posterior of a1,a2,a3 / z1,z2,z3 covers truth;
#   (b) SURFACE level        -- fitted CTmax(x1,x2), z(x1,x2) on an (x1,x2) grid
#                               (engine eval via tls_eval_subpars + custom
#                               newdata) covers the analytical truth at each node;
#   (c) INTERACTION detect.  -- the x1:x2 CrI EXCLUDES 0 when planted non-zero,
#                               and a null companion (a3=z3=0) INCLUDES 0.
# All gated behind RUN_BRMS_TESTS.

# Posterior 95% CrI for a named fixed effect (modeled scale), via brms draws.
# `probs` defaults to the 95% CrI (used by the "excludes 0" signal/power check);
# the truth-in-CrI recovery checks and the null straddles-0 checks pass a WIDE
# interval (0.001/0.999) so a single fitted fixture does not flake across
# platforms (CI's cmdstan draws differ from a locally-cached fixture). Reading
# the two quantile columns by position keeps this robust to brms' probs-based
# column naming.
fixef_ci <- function(fit, par, probs = c(0.025, 0.975)) {
  fx <- brms::fixef(get_brmsfit_or(fit), probs = probs)
  stopifnot(par %in% rownames(fx))
  c(lower = unname(fx[par, 3L]), est = unname(fx[par, "Estimate"]),
    upper = unname(fx[par, 4L]))
}
get_brmsfit_or <- function(x) if (inherits(x, "brmsfit")) x else get_brmsfit(x)

# Per-draw CTmax(x1,x2) and z(x1,x2) at a grid node, via the SAME engine the
# package uses: evaluate the relative mid at temp_c = 0, -h, +h and read
# z = -1/slope, CTmax = T_bar + (log10_tref - mid0)/slope.
surface_at <- function(wf, grid, h = 1e-3) {
  fit <- get_brmsfit(wf); Tbar <- wf$meta$temp_mean
  nd <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
    data.frame(temp_c = c(0, -h, h), x1 = grid$x1[i], x2 = grid$x2[i])))
  for (v in setdiff(names(fit$data), names(nd))) nd[[v]] <- fit$data[[v]][1]
  sp <- tls_eval_subpars(fit, nd, wf$meta$bounds, mode = "relative")
  lapply(seq_len(nrow(grid)), function(i) {
    cols  <- (i - 1) * 3 + 1:3
    mid0  <- sp$mid[, cols[1]]; slope <- (sp$mid[, cols[3]] - sp$mid[, cols[2]]) / (2 * h)
    list(z  = -1 / slope,
         ct = Tbar + (wf$meta$log10_tref - mid0) / slope)
  })
}

test_that("PRE-FLIGHT: fit_4pl accepts ~ x1*x2 on ctmax/z and exposes the interaction coefficients", {
  skip_unless_brms()
  wf <- load_fixture_workflow_cont2()
  expect_equal(wf$meta$parameterization, "direct")
  # fit_4pl records the moderator names (x1, x2) in meta$group_vars regardless of
  # whether they are factor or continuous -- this is what lets the engine cross
  # their levels for per-condition reads. The coefficient-level recovery below
  # works directly off the brms fixed effects, so this does not depend on the
  # group-var bookkeeping.
  expect_true(all(c("x1", "x2") %in% wf$meta$group_vars))
  fx <- rownames(brms::fixef(get_brmsfit(wf)))
  expect_true(all(c("CTmaxdev_Intercept", "CTmaxdev_x1", "CTmaxdev_x2", "CTmaxdev_x1:x2",
                    "logz_Intercept", "logz_x1", "logz_x2", "logz_x1:x2") %in% fx))
})

test_that("(a) coefficient-level: a1,a2,a3 (CTmaxdev) and z1,z2,z3 (logz) recover truth in CrI", {
  skip_unless_brms()
  wf <- load_fixture_workflow_cont2()
  tr <- attr(wf, "truth_cont2")
  fit <- get_brmsfit(wf)

  # Truth-in-CrI on a single fixture: use a WIDE (99.8%) interval (see fixef_ci note).
  wide <- c(0.001, 0.999)
  # CTmaxdev coefficients (deg C) -- planted a0..a3.
  for (nm in c("a0", "a1", "a2", "a3")) {
    par <- paste0("CTmaxdev_", switch(nm, a0 = "Intercept", a1 = "x1", a2 = "x2", a3 = "x1:x2"))
    ci  <- fixef_ci(fit, par, probs = wide)
    expect_gte(tr$a[[nm]], ci["lower"]); expect_lte(tr$a[[nm]], ci["upper"])
  }
  # logz coefficients (log scale) -- planted z0..z3.
  for (nm in c("z0", "z1", "z2", "z3")) {
    par <- paste0("logz_", switch(nm, z0 = "Intercept", z1 = "x1", z2 = "x2", z3 = "x1:x2"))
    ci  <- fixef_ci(fit, par, probs = wide)
    expect_gte(tr$zc[[nm]], ci["lower"]); expect_lte(tr$zc[[nm]], ci["upper"])
  }
})

test_that("(b) surface-level: fitted CTmax(x1,x2) and z(x1,x2) cover the truth surface at each grid node", {
  skip_unless_brms()
  wf <- load_fixture_workflow_cont2()
  tr <- attr(wf, "truth_cont2")
  grid <- expand.grid(x1 = c(-1, 0, 1), x2 = c(-1, 0, 1))
  surf <- surface_at(wf, grid)

  for (i in seq_len(nrow(grid))) {
    x1 <- grid$x1[i]; x2 <- grid$x2[i]
    z_true  <- exp(tr$zc["z0"] + tr$zc["z1"] * x1 + tr$zc["z2"] * x2 + tr$zc["z3"] * x1 * x2)
    ct_true <- tr$T_bar + (tr$a["a0"] + tr$a["a1"] * x1 + tr$a["a2"] * x2 + tr$a["a3"] * x1 * x2)
    # wide (99.8%) interval: single-fixture surface coverage, robust across platforms.
    zq  <- stats::quantile(surf[[i]]$z,  c(0.001, 0.999), names = FALSE, na.rm = TRUE)
    cq  <- stats::quantile(surf[[i]]$ct, c(0.001, 0.999), names = FALSE, na.rm = TRUE)
    expect_gte(unname(z_true),  zq[1]); expect_lte(unname(z_true),  zq[2])
    expect_gte(unname(ct_true), cq[1]); expect_lte(unname(ct_true), cq[2])
  }
  # Sanity: the surface is non-flat -- z at (1,1) differs from z at (-1,-1).
  expect_gt(abs(stats::median(surf[[which(grid$x1 == 1 & grid$x2 == 1)]]$z) -
                stats::median(surf[[which(grid$x1 == -1 & grid$x2 == -1)]]$z)), 0.5)
})

test_that("(c) interaction detectability: x1:x2 CrI EXCLUDES 0 when planted, INCLUDES 0 when null", {
  skip_unless_brms()
  wf  <- load_fixture_workflow_cont2()       # a3 = 0.6, z3 = 0.20 (planted)
  wfn <- load_fixture_workflow_cont2_null()  # a3 = z3 = 0

  ci_ct  <- fixef_ci(get_brmsfit(wf),  "CTmaxdev_x1:x2")
  ci_lz  <- fixef_ci(get_brmsfit(wf),  "logz_x1:x2")
  # Planted interactions are positive and the 95% CrI excludes 0.
  expect_gt(ci_ct["lower"], 0)
  expect_gt(ci_lz["lower"], 0)

  ci_ct0 <- fixef_ci(get_brmsfit(wfn), "CTmaxdev_x1:x2", probs = c(0.001, 0.999))
  ci_lz0 <- fixef_ci(get_brmsfit(wfn), "logz_x1:x2",     probs = c(0.001, 0.999))
  # Null interactions: the wide CrI straddles 0 (robust; the 95% power check above
  # is the strict signal test).
  expect_lt(ci_ct0["lower"], 0); expect_gt(ci_ct0["upper"], 0)
  expect_lt(ci_lz0["lower"], 0); expect_gt(ci_lz0["upper"], 0)
  # And the null's true interaction (0) is inside its CrI.
  expect_gte(0, ci_ct0["lower"]); expect_lte(0, ci_ct0["upper"])
})

test_that("extract_4pl_pars crosses the continuous moderators and yields finite per-temp params", {
  skip_unless_brms()
  wf <- load_fixture_workflow_cont2()
  tr <- attr(wf, "truth_cont2")
  # x1/x2 are recorded as moderators, so the faithful extractor crosses their
  # observed unique values and prepends the moderator columns. It must return
  # finite (low, up, k, mid) at every (covariate-combo x temperature) cell.
  p <- extract_4pl_pars(wf, temps = c(34, 36))
  expect_true(all(c("x1", "x2", ".draw", "temp", "low", "up", "k", "mid") %in% names(p)))
  expect_true(all(is.finite(p$low) & is.finite(p$up) & is.finite(p$k) & is.finite(p$mid)))
  expect_setequal(unique(p$temp), c(34, 36))
  # The asymptotes were planted constant in the moderators -> roughly stable
  # across covariate combos (median up near the planted 0.97).
  expect_lt(abs(stats::median(p$up) - tr$up), 0.05)
})
