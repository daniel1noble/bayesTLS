# Golden fixtures for the model builders: for each configuration in
# helper-golden.R the formulas, the prior table and the brms-generated Stan code
# are frozen under tests/testthat/golden/ (not fixtures/, which .Rbuildignore
# drops from the tarball). Any change to what fit_4pl(fit = FALSE) and its
# siblings produce fails here. This is the "nothing changed" gate for the
# alternative-models work: brms's file_refit = "on_change" compares Stan code,
# so a drifted string would silently refit every cached case-study model
# (scripts/package/check_cached_fits_stancode.R checks those fits directly).
#
# Two layers:
#   * always: the formulas as base R deparses them and the prior table's
#     defining columns (brms-independent);
#   * when the installed brms matches the version recorded in the fixtures:
#     brms's printed formula, the full prior table and the Stan code.
# Skipped on CRAN (and so under R CMD check --as-cran, including CI): it is a
# development gate, run by devtools::test(). Regenerate deliberately with
#   REGENERATE_GOLDEN=true devtools::test(filter = "golden")

test_that("builder outputs match the golden fixtures", {
  skip_on_cran()
  skip_if_not_installed("brms")
  regen <- identical(Sys.getenv("REGENERATE_GOLDEN"), "true")
  cfgs  <- golden_configs()
  for (nm in names(cfgs)) {
    path <- file.path(golden_dir(), paste0(nm, ".rds"))
    now  <- golden_snapshot(cfgs[[nm]]$build())
    if (regen) {
      dir.create(golden_dir(), showWarnings = FALSE, recursive = TRUE)
      saveRDS(now, path)
      succeed()
      next
    }
    expect_true(file.exists(path), info = paste("missing golden fixture", nm))
    if (!file.exists(path)) next
    gold <- readRDS(path)
    expect_identical(now$formula_parts, gold$formula_parts,
                     info = paste(nm, "formulas"))
    expect_identical(now$prior_core, gold$prior_core,
                     info = paste(nm, "prior table"))
    if (identical(now$brms, gold$brms)) {
      expect_identical(now$formula, gold$formula,
                       info = paste(nm, "printed formula"))
      expect_identical(now$prior, gold$prior, info = paste(nm, "full prior"))
      expect_identical(now$stancode, gold$stancode,
                       info = paste(nm, "Stan code"))
    }
  }
})

test_that("the golden configurations cover every model", {
  models <- vapply(golden_configs(), `[[`, "", "model")
  expect_setequal(unique(models), names(tdt_models()))
})
