# NO-DRIFT GATE for the posterior_linpred unification. Re-runs the migrated
# extraction functions over the same fixtures + real case-study fits captured
# BEFORE the refactor (tests/testthat/fixtures/golden_preunify.rds, built by the
# Step-0 snapshot) and asserts the results have NOT moved. Because each function
# kept its exact derive formula and changed only the data source (coefficient
# parsing -> posterior_linpred), every quantity should match to ~machine
# precision. Gated (needs the cached brms fits). If golden_preunify.rds is
# absent, the gate skips with a clear message (rebuild via the Step-0 snapshot).

parity_workflows <- function() {
  list(midpoint = load_fixture_workflow(),
       beta     = load_fixture_workflow_beta(),
       direct   = load_fixture_workflow_direct())
}

test_that("unification did not move tdt_parameter_table / derive_z / predict_survival_curves", {
  skip_unless_brms()
  gp <- here::here("tests", "testthat", "fixtures", "golden_preunify.rds")
  if (!file.exists(gp)) skip("local-only no-drift gate: golden_preunify.rds is gitignored and present only where the Step-0 snapshot was built (it pins the exact cached fits); skipped in CI / fresh checkouts.")
  g <- readRDS(gp); wfs <- parity_workflows()

  for (nm in names(wfs)) {
    expect_lt(max(abs(tdt_parameter_table(wfs[[nm]])$median - g[[nm]]$par_table$median)), 1e-4)
    expect_lt(abs(derive_z(wfs[[nm]], seed = 1)$summary$z_median - g[[nm]]$derive_z$summary$z_median), 1e-4)
    # predict_heat_injury parity is intentionally NOT checked against this golden:
    # its dose integral was corrected after the snapshot was captured (per-interval
    # dt + no endpoint over-count), so it no longer matches the pre-fix values.
    # HI is covered by test-predict_heat_injury.R.
    new_sc <- predict_survival_curves(wfs[[nm]], temps = head(sort(unique(wfs[[nm]]$data$temp)), 3),
                                      durations = c(0.5, 2, 8), ndraws = NULL)$summary
    expect_lt(max(abs(new_sc$survival_median - g[[nm]]$surv$survival_median)), 1e-6)
  }
})
