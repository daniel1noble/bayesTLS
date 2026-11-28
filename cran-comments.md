# cran-comments.md

## Submission

This is an update of bayesTLS from 1.0.0 (on CRAN since 2026-07-21) to 1.3.0.
Versions 1.1.0 and 1.2.0 were released on GitHub only. The update adds
log-normal and Weibull dose-response models (`fit_lognormal()`,
`fit_weibull()`), covariance-matrix (e.g. phylogenetic) random effects, fixed
asymptotes and censored continuous-proportion readings. It also fixes several
bugs in how thresholds are derived from the posterior. NEWS.md lists every
change. No exported function was removed, and every new argument comes last in
its function's signature.

## Test environments

* local: macOS 27.2 (aarch64-apple-darwin23), R 4.6.2. `R CMD check --as-cran`
  on the built tarball, with the PDF manual.
* GitHub Actions: `R CMD check` on ubuntu-latest, R-release, passes on the
  development version this release is built from. The release itself changes
  only DESCRIPTION, NEWS.md and the help pages (two text corrections; the
  rest were regenerated with roxygen2 8.1.1).

## R CMD check results

Local: 0 errors | 0 warnings | 2 notes. Neither is a package issue.

* "Skipping checking HTML validation: 'tidy' doesn't look like recent enough
  HTML Tidy" comes from the old HTML Tidy that ships with macOS. With V8
  installed, the math-rendering check of the HTML manual passes.
* "Examples with CPU (user + system) or elapsed time > 5s" lists
  `make_4pl_formula` at 1.4 s of CPU but 11 s elapsed, because other jobs were
  loading the machine during the check. The same examples raised no timing
  note in an earlier check on an idle machine.

The CRAN incoming feasibility check reports one item, as it did for 1.0.0.

* Suggests or Enhances not in mainstream repositories: cmdstanr. `cmdstanr` is
  available from the Stan development repository, which is declared in
  `Additional_repositories:` (<https://stan-dev.r-universe.dev>), and the check
  confirms "cmdstanr yes". bayesTLS never calls cmdstanr functions itself.
  cmdstanr is one of the two Stan backends that brms can use, and the fitting
  functions also work with `backend = "rstan"`, the RStan backend that brms
  imports.

If the check lists possibly misspelled words in DESCRIPTION (CTmax, crit,
Schoolfield, Nakagawa), none are misspellings. "CTmax" is the standard
abbreviation for the critical thermal maximum, "crit" is part of "T_crit",
"Schoolfield" refers to the Sharpe-Schoolfield rate model and "Nakagawa" is a
co-author's surname.

## Reverse dependencies

There are no reverse dependencies on CRAN
(`tools::package_dependencies("bayesTLS", reverse = TRUE, which = "all")`,
checked on 2026-11-28).

## Notes for the reviewer (unchanged since 1.0.0)

* Examples that fit models use `\dontrun{}`. Fitting even a small model needs
  Stan, a C++ toolchain and far longer than the few seconds CRAN allows for an
  example, so these examples cannot run during the check. All examples that do
  not need a fitted model run normally.
* The integration tests that fit Stan models are gated behind the
  `RUN_BRMS_TESTS` environment variable and skip by default, so the check runs
  quickly and does not invoke a compiler. The cached fits those tests use are
  excluded from the tarball through `.Rbuildignore`.
