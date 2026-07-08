# PILLAR 2 -- two-group (categorical) recovery across families.
#
# These EXTEND the existing 2-group tests (test-direct-grouped.R and
# test-midpoint-grouped.R, which use the beta-binomial simulate_tdt_2group()
# data) to the two families those do not cover -- PURE BINOMIAL and continuous
# BETA -- and add explicit between-group CONTRAST (delta-z, delta-CTmax)
# truth-in-CrI checks and a direct==midpoint per-group agreement check on the
# same binomial data. All gated behind RUN_BRMS_TESTS.
#
# Each group has KNOWN, deliberately-distinct (z, CTmax): a large gap means a
# swapped group->parameter mapping cannot pass by coincidence.

# Per-draw per-group contrast tibble for a derived quantity (z or CTmax draws),
# joined on .draw so the contrast is genuinely paired within draw.
group_contrast <- function(draws, value_col, gcol = "grp") {
  a <- draws[draws[[gcol]] == "A", c(".draw", value_col)]
  b <- draws[draws[[gcol]] == "B", c(".draw", value_col)]
  names(a)[2] <- "vA"; names(b)[2] <- "vB"
  m <- merge(a, b, by = ".draw")
  m$delta <- m$vB - m$vA
  m
}

# Shared assertions: per-group z & CTmax recovery, correct mapping (large gap
# recovered with the right sign), and contrast (delta) recovery.
#
# Truth-IN-INTERVAL checks use a WIDE (99.8%) credible interval, NOT the 95% CrI.
# On a single fitted fixture, 95%-CrI coverage of a planted truth is a ~5%
# coin-flip per bound, and with ~16 such bounds it flakes across platforms (CI's
# cmdstan build produces slightly different draws than the locally-cached
# fixture). Recovery ACCURACY is pinned by the point-estimate tolerance and
# mapping checks below; calibrated interval coverage is validated by the
# manuscript's coverage simulation, not a single-fixture unit test. The
# "excludes 0" distinctness check keeps the standard 95% interval (it is a
# signal test, not a coverage test).
assert_two_group_recovery <- function(wf, tr, z_tol = 2.0, ct_tol = 1.2) {
  ci <- c(0.001, 0.999)
  out <- extract_tdt(wf, t_ref = 60, ndraws = NULL)
  zs <- out$z$summary; cs <- out$CTmax$summary
  zd <- out$z$draws;   cd <- out$CTmax$draws
  expect_true("grp" %in% names(zs))
  expect_setequal(as.character(zs$grp), c("A", "B"))

  for (g in c("A", "B")) {
    zg <- zs[zs$grp == g, ]; cg <- cs[cs$grp == g, ]
    # per-group point estimate near truth ...
    expect_lt(abs(zg$z_median        - tr[[g]]$z),         z_tol)
    expect_lt(abs(cg$temp_median     - tr[[g]]$CTmax_1hr), ct_tol)
    # ... and truth inside a wide (99.8%) credible interval (see note above).
    qzg <- stats::quantile(zd$z[zd$grp == g],    ci, names = FALSE)
    qcg <- stats::quantile(cd$temp[cd$grp == g], ci, names = FALSE)
    expect_gte(tr[[g]]$z,         qzg[1]); expect_lte(tr[[g]]$z,         qzg[2])
    expect_gte(tr[[g]]$CTmax_1hr, qcg[1]); expect_lte(tr[[g]]$CTmax_1hr, qcg[2])
  }

  # Mapping guard: A has the SHALLOWER slope (larger z) and HIGHER CTmax by
  # construction, so z_A > z_B and CTmax_A > CTmax_B must hold in the posterior.
  zA <- zs$z_median[zs$grp == "A"]; zB <- zs$z_median[zs$grp == "B"]
  cA <- cs$temp_median[cs$grp == "A"]; cB <- cs$temp_median[cs$grp == "B"]
  expect_gt(zA - zB, 0.8)
  expect_gt(cA - cB, 0.2)

  # Between-group CONTRASTS, paired within draw: truth inside the wide CrI.
  dz <- group_contrast(out$z$draws, "z")
  dc <- group_contrast(out$CTmax$draws, "temp")
  truth_dz <- tr$B$z - tr$A$z
  truth_dc <- tr$B$CTmax_1hr - tr$A$CTmax_1hr
  qz <- stats::quantile(dz$delta, ci, names = FALSE)
  qc <- stats::quantile(dc$delta, ci, names = FALSE)
  expect_gte(truth_dz, qz[1]); expect_lte(truth_dz, qz[2])
  expect_gte(truth_dc, qc[1]); expect_lte(truth_dc, qc[2])
  # Delta-z 95% CrI excludes 0 (groups genuinely distinct) -- a signal test, so
  # it keeps the standard 95% interval.
  expect_lt(stats::quantile(dz$delta, 0.975, names = FALSE), 0)
  invisible(out)
}

test_that("PURE BINOMIAL 2-group (direct): per-group z/CTmax, mapping, contrasts recovered", {
  skip_unless_brms()
  wf <- load_fixture_workflow_2group_binom_direct()
  expect_equal(wf$meta$parameterization, "direct")
  tr <- attr(wf, "truth_2group")
  assert_two_group_recovery(wf, tr)
})

test_that("PURE BINOMIAL 2-group: per-group upper asymptote recovered (A's > B's)", {
  skip_unless_brms()
  wf <- load_fixture_workflow_2group_binom_direct()
  tr <- attr(wf, "truth_2group")
  tab <- tdt_parameter_table(wf, by = "grp")     # per-group low/up/k/CTmax/z
  upA <- tab$median[tab$grp == "A" & grepl("^up", tab$parameter)]
  upB <- tab$median[tab$grp == "B" & grepl("^up", tab$parameter)]
  # Truth: up_A = 0.97, up_B = 0.90 -> A's plateau is higher and each near truth.
  expect_lt(abs(upA - tr$A$up), 0.08)
  expect_lt(abs(upB - tr$B$up), 0.08)
  expect_gt(upA - upB, 0.0)
})

test_that("BETA continuous 2-group (direct): per-group z/CTmax, mapping, contrasts recovered", {
  skip_unless_brms()
  wf <- load_fixture_workflow_2group_beta_direct()
  expect_equal(wf$meta$parameterization, "direct")
  tr <- attr(wf, "truth_2group")
  assert_two_group_recovery(wf, tr)
})

test_that("BINOMIAL 2-group: direct (~0+grp) and midpoint (by=grp) agree per group", {
  skip_unless_brms()
  wf_d <- load_fixture_workflow_2group_binom_direct()
  wf_m <- load_fixture_workflow_2group_binom_mid()
  expect_equal(wf_m$meta$parameterization, "midpoint")
  ed <- extract_tdt(wf_d, t_ref = 60, ndraws = NULL)
  em <- extract_tdt(wf_m, t_ref = 60, ndraws = NULL)
  for (g in c("A", "B")) {
    zd <- ed$z$summary$z_median[ed$z$summary$grp == g]
    zm <- em$z$summary$z_median[em$z$summary$grp == g]
    cd <- ed$CTmax$summary$temp_median[ed$CTmax$summary$grp == g]
    cm <- em$CTmax$summary$temp_median[em$CTmax$summary$grp == g]
    # Same data, two parameterisations -> medians coincide up to MCMC noise.
    expect_lt(abs(zd - zm), 0.8, label = paste("z group", g))
    expect_lt(abs(cd - cm), 0.5, label = paste("CTmax group", g))
  }
})

test_that("BINOMIAL 2-group: tls(by=) and extract_tdt(by=) agree per group", {
  skip_unless_brms()
  wf <- load_fixture_workflow_2group_binom_direct()
  et <- extract_tdt(wf, t_ref = 60, ndraws = NULL)
  tl <- tls(wf, by = "grp", params = c("z", "ctmax"), mode = "relative",
            t_ref = 60, temp_mean = wf$meta$temp_mean)$summary
  for (g in c("A", "B")) {
    z_et <- et$z$summary$z_median[et$z$summary$grp == g]
    z_tl <- tl$median[tl$grp == g & tl$quantity == "z"]
    expect_lt(abs(z_et - z_tl), 0.3)
  }
})
