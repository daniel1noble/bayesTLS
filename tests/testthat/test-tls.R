# Accuracy tests for tls(). The brms-free block checks argument validation; the
# gated block validates numerical correctness against the exact -1/beta1
# identity, the simulation truth, and extract_tdt(), plus the per-group machinery.

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

test_that("tls() matches extract_tdt() for z and CTmax (relative and absolute)", {
  skip_unless_brms()
  wf     <- load_fixture_workflow()
  nd_all <- brms::ndraws(get_brmsfit(wf))
  g <- function(s, q) s$median[s$quantity == q]

  for (m in c("relative", "absolute")) {
    et <- extract_tdt(wf, target_surv = m, t_ref = 60, ndraws = nd_all)
    tl <- tls(wf, params = c("z", "ctmax"), mode = m, t_ref = 60)$summary
    expect_equal(g(tl, "z"),     et$z$summary$z_median,        tolerance = 0.05,
                 info = paste("z,", m))
    expect_equal(g(tl, "CTmax"), et$CTmax$summary$temp_median, tolerance = 0.1,
                 info = paste("CTmax,", m))
  }
})

test_that("tls() and extract_tdt() share one engine: identical z/CTmax on a common grid", {
  skip_unless_brms()
  wf     <- load_fixture_workflow()
  nd_all <- brms::ndraws(get_brmsfit(wf))
  tbar   <- wf$meta$temp_mean
  obs    <- sort(unique(wf$data$temp))          # extract_tdt pools z over these
  gz  <- function(o) o$summary$median[o$summary$quantity == "z"]
  gc  <- function(o) o$summary$median[o$summary$quantity == "CTmax"]
  fin <- function(o, q) { d <- o$draws[o$draws$quantity == q, c(".draw", "value")]
                          d[is.finite(d$value), ] }

  # Common grid = observed assay temps for BOTH z-pooling and CTmax inversion, so
  # a single tls() call matches a single extract_tdt() call for z AND CTmax. A
  # revert to the old global least-squares slope would diverge on the (bent)
  # absolute curve far beyond 1e-6.
  for (m in c("relative", "absolute")) {
    et <- extract_tdt(wf, target_surv = m, t_ref = 60, temp_grid = obs, ndraws = nd_all)
    tl <- tls(wf, params = c("z", "ctmax"), target_surv = m, t_ref = 60,
              temp_grid = obs - tbar)
    expect_equal(gz(tl), et$z$summary$z_median,        tolerance = 1e-6,
                 info = paste("z median,", m))
    expect_equal(gc(tl), et$CTmax$summary$temp_median, tolerance = 1e-6,
                 info = paste("CTmax median,", m))
    # per-draw parity (all draws, in order -> aligned 1:1)
    zt <- merge(fin(tl, "z"),     et$z$draws[c(".draw", "z")],        by = ".draw")
    ct <- merge(fin(tl, "CTmax"), et$CTmax$draws[c(".draw", "temp")], by = ".draw")
    expect_lt(max(abs(zt$value - zt$z)),    1e-6)
    expect_lt(max(abs(ct$value - ct$temp)), 1e-6)
  }
})

test_that("tls(absolute) == extract_tdt(absolute) on a strongly bent curve (LS-slope-fix guard)", {
  skip_unless_brms()
  wf   <- load_bent_workflow()                  # skips if the cached bent fit is absent
  tbar <- wf$meta$temp_mean
  obs  <- sort(unique(wf$data$temp))
  et <- extract_tdt(wf, target_surv = "absolute", t_ref = 60, time_multiplier = 1,
                    temp_grid = obs, ndraws = NULL)
  tl <- tls(wf, params = c("z", "ctmax"), target_surv = "absolute", t_ref = 60,
            time_multiplier = 1, temp_grid = obs - tbar)
  gz <- tl$summary$median[tl$summary$quantity == "z"]
  gc <- tl$summary$median[tl$summary$quantity == "CTmax"]
  expect_equal(gz, et$z$summary$z_median,        tolerance = 1e-6)   # bent curve, still identical
  expect_equal(gc, et$CTmax$summary$temp_median, tolerance = 1e-6)
  zt <- merge(tl$draws[tl$draws$quantity == "z",     c(".draw", "value")],
              et$z$draws[c(".draw", "z")], by = ".draw")
  ct <- merge(tl$draws[tl$draws$quantity == "CTmax", c(".draw", "value")],
              et$CTmax$draws[c(".draw", "temp")], by = ".draw")
  expect_lt(max(abs(zt$value - zt$z),    na.rm = TRUE), 1e-6)
  expect_lt(max(abs(ct$value - ct$temp), na.rm = TRUE), 1e-6)
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
