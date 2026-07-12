# cran-comments.md

## Submission

This is a new submission: the first release of **bayesTLS** (1.0.0).

## Test environments

* local: macOS 26.5.2 (aarch64-apple-darwin20), R 4.4.2 — `R CMD check --as-cran`

<!--
BEFORE SUBMITTING, run and record the results of at least:

  * win-builder, both release and devel:
        devtools::check_win_release()
        devtools::check_win_devel()
  * R-hub (Linux + Windows):
        rhub::rhub_check()

Add each platform and its result below. Do not submit until these have
actually been run — CRAN will run them regardless.
-->

## R CMD check results

Local check: **0 errors | 0 warnings | 1 note**.

The one NOTE is from `checking CRAN incoming feasibility`:

```
The Title field should be in title case. Current version is:
  'Joint Bayesian 4PL Models for Thermal Load Sensitivity'
In title case that is:
  'Joint Bayesian 4pl Models for Thermal Load Sensitivity'
```

This is a false positive. "4PL" is the standard acronym for the
*four-parameter logistic* model and is correctly capitalised; lowercasing it to
"4pl" would be wrong. The Title is otherwise in title case.

## Notes for the reviewer

* **Examples that fit models use `\dontrun{}`.** `bayesTLS` fits its models with
  Stan (via `brms`). Fitting even a deliberately small model requires a working
  C++ toolchain and takes far longer than the few seconds CRAN allots to an
  example, so the model-fitting examples cannot be run during the check. All
  examples that do *not* require a fit run normally and are checked. This
  follows the convention used by other Stan-backed packages.

* **`cmdstanr` is in `Suggests` and is not on CRAN.** It is available from the
  Stan development repository, which is declared in `Additional_repositories:`
  (<https://stan-dev.r-universe.dev>). The package works without it — `brms`
  falls back to the RStan backend — and all code paths that touch `cmdstanr` are
  guarded with `requireNamespace()`.

* **Cached model fits are not shipped.** The integration tests that exercise real
  Stan fits are gated behind the `RUN_BRMS_TESTS` environment variable and skip
  by default, so the check runs quickly and does not invoke a compiler. The
  cached fits those tests rely on are excluded from the tarball via
  `.Rbuildignore`.
