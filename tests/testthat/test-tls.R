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
