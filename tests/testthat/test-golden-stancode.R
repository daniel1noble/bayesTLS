# Golden fixtures for the model builders. For each canonical configuration the
# formula (as printed), the prior table and the brms-generated Stan code are
# frozen under fixtures/models/. Any change to what make_4pl_formula() /
# make_4pl_priors() / fit_4pl(fit = FALSE) produce fails here. This is the
# "nothing changed for the log-logistic" gate for the alternative-models work:
# brms's file_refit = "on_change" compares Stan code, so a drifted string would
# silently refit every cached case-study model.
#
# Regenerate deliberately (after an intended change) with
#   REGENERATE_GOLDEN=true devtools::test(filter = "golden")
# Stan code is compared only when the brms version matches the one recorded in
# the fixture; the formula and prior table are compared always.

golden_dir <- testthat::test_path("fixtures", "models")

g_std <- function() {
  standardize_data(simulate_tdt(seed = 1), temp = "T", duration = "t_hours",
                   n_total = "n_trials", n_surv = "y_alive",
                   duration_unit = "hours")
}
g_std2 <- function() {
  standardize_data(simulate_tdt_2group(seed = 2), temp = "T",
                   duration = "t_hours", n_total = "n_trials",
                   n_surv = "y_alive", duration_unit = "hours")
}
g_beta <- function() {
  standardize_data(simulate_tdt_beta(seed = 3), temp = "T",
                   duration = "t_hours", proportion = "y_prop",
                   duration_unit = "hours")
}

golden_configs <- list(
  mid_bb_default = function() fit_4pl(g_std(), fit = FALSE),
  mid_bb_by      = function() fit_4pl(g_std2(), by = "grp", fit = FALSE),
  mid_bb_shape   = function() fit_4pl(g_std(), low = ~ 1, up = ~ temp_c,
                                      k = ~ temp_c, mid = ~ temp_c, fit = FALSE),
  mid_binomial   = function() fit_4pl(g_std(),
                                      family = stats::binomial(link = "identity"),
                                      fit = FALSE),
  mid_beta       = function() fit_4pl(g_beta(), fit = FALSE),
  direct_rel_1   = function() fit_4pl(g_std(), ctmax = ~ 1, z = ~ 1, fit = FALSE),
  direct_rel_grp = function() fit_4pl(g_std2(), ctmax = ~ 0 + grp, z = ~ 0 + grp,
                                      fit = FALSE),
  direct_abs_1   = function() fit_4pl(g_std(), ctmax = ~ 1, z = ~ 1,
                                      threshold = "absolute", fit = FALSE),
  mid_fixed_asy  = function() fit_4pl(g_std(), fix_asymptotes = TRUE, fit = FALSE)
)

snapshot_of <- function(wf) {
  list(formula  = utils::capture.output(print(wf$formula)),
       prior    = as.data.frame(wf$prior),
       # suppressWarnings: brms reports unused global `b` priors for
       # intercept-only direct configs (pre-existing package behaviour); the
       # Stan code string is what is frozen here.
       stancode = suppressWarnings(brms::make_stancode(wf$formula, data = wf$data,
                                                       prior = wf$prior)),
       brms     = as.character(utils::packageVersion("brms")))
}

test_that("builder outputs match the golden fixtures", {
  skip_if_not_installed("brms")
  regen <- identical(Sys.getenv("REGENERATE_GOLDEN"), "true")
  for (nm in names(golden_configs)) {
    path <- file.path(golden_dir, paste0(nm, ".rds"))
    now  <- snapshot_of(golden_configs[[nm]]())
    if (regen) {
      dir.create(golden_dir, showWarnings = FALSE, recursive = TRUE)
      saveRDS(now, path)
      succeed()
      next
    }
    expect_true(file.exists(path), info = paste("missing golden fixture", nm))
    if (!file.exists(path)) next
    gold <- readRDS(path)
    expect_identical(now$formula, gold$formula, info = paste(nm, "formula"))
    expect_identical(now$prior,   gold$prior,   info = paste(nm, "prior"))
    if (identical(now$brms, gold$brms)) {
      expect_identical(now$stancode, gold$stancode,
                       info = paste(nm, "Stan code"))
    } else {
      testthat::skip(sprintf(
        "brms %s differs from the golden fixtures' %s: Stan code not compared",
        now$brms, gold$brms))
    }
  }
})
