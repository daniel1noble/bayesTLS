# cran-comments.md

## Submission

This is a new submission: the first release of **bayesTLS** (1.0.0).

## Test environments

* local: macOS (aarch64-apple-darwin), R 4.5 — `R CMD check --as-cran` (with
  PDF manual)
* win-builder, R-release (R 4.6.1) — via `devtools::check_win_release()`
* win-builder, R-devel — via `devtools::check_win_devel()`
* R-hub, Linux (R-devel) container — via `rhub::rhub_check()`
* R-hub, Windows (R-devel) — via `rhub::rhub_check()`

## R CMD check results

All environments: **0 errors | 0 warnings**.

win-builder (release and devel) each report **1 NOTE**; R-hub Linux and Windows
pass. The one NOTE is `checking CRAN incoming feasibility`, and every item in it
is expected for a first submission:

* **New submission.** This is the first release.

* **Possibly misspelled words in DESCRIPTION: CTmax, crit, Schoolfield,
  Nakagawa.** None are misspellings: "CTmax" is the standard critical-thermal-
  maximum abbreviation, "crit" is part of "T_crit", "Schoolfield" refers to the
  Sharpe-Schoolfield rate model, and "Nakagawa" is a co-author's surname.

* **Suggests or Enhances not in mainstream repositories: cmdstanr.** `cmdstanr`
  is available from the Stan development repository, declared in
  `Additional_repositories:` (<https://stan-dev.r-universe.dev>); the check
  confirms "cmdstanr yes". The package works without it — `brms` falls back to
  the RStan backend — and all code paths that touch `cmdstanr` are guarded with
  `requireNamespace()`.

* **Title field should be in title case ("4PL" -> "4pl").** A false positive.
  "4PL" is the standard acronym for the *four-parameter logistic* model and is
  correctly capitalised; lowercasing it would be wrong.

## Notes for the reviewer

* **Examples that fit models use `\dontrun{}`.** `bayesTLS` fits its models with
  Stan (via `brms`). Fitting even a deliberately small model requires a working
  C++ toolchain and takes far longer than the few seconds CRAN allots to an
  example, so the model-fitting examples cannot be run during the check. All
  examples that do *not* require a fit run normally and are checked. This
  follows the convention used by other Stan-backed packages.

* **Cached model fits are not shipped.** The integration tests that exercise real
  Stan fits are gated behind the `RUN_BRMS_TESTS` environment variable and skip
  by default, so the check runs quickly and does not invoke a compiler. The
  cached fits those tests rely on are excluded from the tarball via
  `.Rbuildignore`.
