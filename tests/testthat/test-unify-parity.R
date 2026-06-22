# NO-DRIFT GATE for the posterior_linpred unification. Re-runs the migrated
# extraction functions over the same fixtures + real case-study fits captured
# BEFORE the refactor (tests/testthat/fixtures/golden_preunify.rds, built by the
# Step-0 snapshot) and asserts the results have NOT moved. Because each function
# kept its exact derive formula and changed only the data source (coefficient
# parsing -> posterior_linpred), every quantity should match to ~machine
# precision. Gated (needs the cached brms fits). If golden_preunify.rds is
# absent, the gate skips with a clear message (rebuild via the Step-0 snapshot).

cb <- function(...) getFromNamespace("compute_4pl_bounds", "bayesTLS")(...)

# Wrap a cached brmsfit into a bayes_tls workflow (verifying temp_c alignment),
# identical to the Step-0 snapshot construction.
parity_wrap <- function(fit, std, du = "hours", group_vars = character(0)) {
  meta <- attr(std, "tdt_meta"); tm <- tdt_unit_to_minutes(du)
  stopifnot(isTRUE(all.equal(sort(unique(round(std$temp_c, 6))),
                             sort(unique(round(fit$data$temp_c, 6))))))
  structure(list(fit = fit, data = std,
    meta = utils::modifyList(meta, list(
      bounds = cb(0, 1), lower = 0, upper = 1, parameterization = "midpoint",
      threshold = "relative", t_ref = 60, log10_tref = log10(60 / tm),
      duration_unit = du, group_vars = group_vars,
      grouped = length(group_vars) > 0))), class = "bayes_tls")
}

parity_workflows <- function() {
  wfs <- list(midpoint = load_fixture_workflow(),
              beta     = load_fixture_workflow_beta(),
              direct   = load_fixture_workflow_direct())
  shp <- file.path(here::here("output", "models"), "fit_shrimp_lethal_4pl.rds")
  if (file.exists(shp)) {
    ss <- standardize_data(shrimp_lethal, temp = "Temperature_assay",
            duration = "Duration_exposure_hours", n_total = "N_individuals_after_trial",
            mortality = "Mortality_after_trial", duration_unit = "hours")
    wfs$shrimp <- parity_wrap(readRDS(shp), ss)
  }
  wfs
}

test_that("unification did not move z / CTmax / T_crit (extract_tdt) on any fixture or case study", {
  skip_unless_brms()
  gp <- here::here("tests", "testthat", "fixtures", "golden_preunify.rds")
  if (!file.exists(gp)) skip("local-only no-drift gate: golden_preunify.rds is gitignored and present only where the Step-0 snapshot was built (it pins the exact cached fits); skipped in CI / fresh checkouts.")
  g <- readRDS(gp); wfs <- parity_workflows()

  for (nm in names(wfs)) {
    old <- g[[nm]]$extract
    leth <- !is.null(old$T_crit)
    new <- suppressMessages(extract_tdt(wfs[[nm]], ndraws = 1000, lethal = leth, seed = 1))
    expect_lt(abs(new$z$summary$z_median        - old$z$summary$z_median),        1e-4)
    expect_lt(abs(new$CTmax$summary$temp_median - old$CTmax$summary$temp_median), 1e-4)
    if (leth)
      expect_lt(abs(new$T_crit$summary$temp_median - old$T_crit$summary$temp_median), 0.05)
  }
})

test_that("unification did not move tdt_parameter_table / derive_z / predict_heat_injury / predict_survival_curves", {
  skip_unless_brms()
  gp <- here::here("tests", "testthat", "fixtures", "golden_preunify.rds")
  if (!file.exists(gp)) skip("local-only no-drift gate: golden_preunify.rds is gitignored and present only where the Step-0 snapshot was built (it pins the exact cached fits); skipped in CI / fresh checkouts.")
  g <- readRDS(gp); wfs <- parity_workflows()

  for (nm in names(wfs)) {
    expect_lt(max(abs(tdt_parameter_table(wfs[[nm]])$median - g[[nm]]$par_table$median)), 1e-4)
    expect_lt(abs(derive_z(wfs[[nm]], seed = 1)$summary$z_median - g[[nm]]$derive_z$summary$z_median), 1e-4)
    tr <- data.frame(time = 0:6, temp = seq(min(wfs[[nm]]$data$temp), max(wfs[[nm]]$data$temp), length.out = 7))
    new_hi <- suppressMessages(predict_heat_injury(tr, wfs[[nm]], ndraws = 1000, seed = 1))$summary
    expect_lt(max(abs(new_hi$surv_median - g[[nm]]$hi$surv_median)), 1e-4)
    new_sc <- predict_survival_curves(wfs[[nm]], temps = head(sort(unique(wfs[[nm]]$data$temp)), 3),
                                      durations = c(0.5, 2, 8), ndraws = NULL)$summary
    expect_lt(max(abs(new_sc$survival_median - g[[nm]]$surv$survival_median)), 1e-6)
  }
})

test_that("unification did not move the grouped tls() case study (zebrafish per-group z)", {
  skip_unless_brms()
  gp <- here::here("tests", "testthat", "fixtures", "golden_preunify.rds")
  zfp <- file.path(here::here("output", "models"), "fit_zf_joint_4pl.rds")
  if (!file.exists(gp) || !file.exists(zfp)) skip("local-only no-drift gate: golden snapshot or cached zf fit absent (gitignored); skipped in CI.")
  g <- readRDS(gp)
  if (!is.data.frame(g$zf_tls)) skip("golden zf_tls not captured.")
  zs <- standardize_data(zebrafish_lethal, temp = "assay_temp", duration = "duration_h",
          n_total = "n_total", n_surv = "n_surv", duration_unit = "hours")
  zwf <- parity_wrap(readRDS(zfp), zs, du = "hours", group_vars = "life_stage")
  new <- tls(zwf, by = "life_stage", lethal = TRUE, seed = 1)$summary
  m <- merge(new, g$zf_tls, by = c("life_stage", "quantity"), suffixes = c(".new", ".old"))
  expect_lt(max(abs(m$median.new - m$median.old)), 0.05)   # tls() unchanged by the engine factor-out
})
