# Gated integration tests for covariance-matrix random effects (issue #12):
# a real gr(species, cov = A) fit against committed truth, per-species
# extraction through the tls() re_formula auto-rule, and diagonal parity
# against a plain (1 | species) fit. Fast-tier counterparts (parsing,
# validation, auto-rule structure) live in test-gr-formulas.R.
# Loaders + committed truth: helper-skip.R / fixtures/phylo_truth.rds.

test_that("phylo fit builds, samples, and recovers truth; tls(by) auto-includes
           the phylo random effect (T9 + T10)", {
  skip_unless_brms()
  wf    <- load_fixture_workflow_phylo()
  truth <- attr(wf, "phylo_truth")

  # The fitted brms model carries the covariance term.
  expect_match(paste(deparse(wf$fit$formula$pforms$CTmaxdev), collapse = " "),
               "gr(species, cov = A)", fixed = TRUE)

  # Population-level recovery (single fit: loose gates, no calibration claim).
  glob <- tls(wf, re_formula = NA)$summary
  z_row <- glob[glob$quantity == "z", ]
  expect_true(z_row$lower <= truth$z && truth$z <= z_row$upper)

  # by = "species" is RE-only here -> the auto-rule must announce inclusion.
  # temp_grid extends past the assay range (37-43 C): sp5's true CTmax is
  # ~44.8 C and the CTmax inversion needs the grid to bracket the crossing
  # (otherwise tls() returns honest NAs with a warning -- by design).
  expect_message(
    res <- tls(wf, by = "species", temp_grid = seq(-5, 6, by = 0.25)),
    "random effect")

  s_ct <- res$summary[res$summary$quantity == "CTmax", ]
  expect_identical(nrow(s_ct), length(truth$species))

  # Per-species recovery: BLUP medians track the committed true CTmax values.
  true_ct <- truth$ctmax_global + truth$ctmax_dev
  est     <- s_ct$median[match(names(true_ct), s_ct$species)]
  expect_gt(cor(est, true_ct, method = "spearman"), 0.6)
  expect_lt(max(abs(est - true_ct)), 2.5)   # shrinkage-tolerant absolute gate

  # And the species must genuinely DIFFER (the silent-collapse trap this
  # feature exists to prevent: with re_formula = NA all medians coincide).
  expect_gt(diff(range(est)), 1)
  s_na <- tls(wf, by = "species", re_formula = NA,
              temp_grid = seq(-5, 6, by = 0.25))$summary
  s_na <- s_na[s_na$quantity == "CTmax", ]
  expect_lt(diff(range(s_na$median)), 1e-8) # explicit NA still collapses (override honoured)
})

test_that("gr(species, cov = I) agrees with plain (1 | species) (T11)", {
  skip_unless_brms()
  wf_i <- load_fixture_workflow_phylo_ident()
  wf_p <- load_fixture_workflow_phylo_plain()

  ct_i <- suppressMessages(tls(wf_i, by = "species",
                              temp_grid = seq(-5, 6, by = 0.25)))$summary
  ct_p <- suppressMessages(tls(wf_p, by = "species",
                              temp_grid = seq(-5, 6, by = 0.25)))$summary
  ct_i <- ct_i[ct_i$quantity == "CTmax", ]
  ct_p <- ct_p[ct_p$quantity == "CTmax", ]
  m <- merge(ct_i, ct_p, by = "species", suffixes = c(".ident", ".plain"))

  # Same model up to the (identity) covariance injection: agreement within
  # MCMC noise (different Stan programs, so not bit-identical).
  expect_gt(cor(m$median.ident, m$median.plain), 0.95)
  expect_lt(max(abs(m$median.ident - m$median.plain)), 0.6)

  z_i <- suppressMessages(tls(wf_i))$summary
  z_p <- suppressMessages(tls(wf_p))$summary
  expect_lt(abs(z_i$median[z_i$quantity == "z"] -
                z_p$median[z_p$quantity == "z"]), 0.25)
})
