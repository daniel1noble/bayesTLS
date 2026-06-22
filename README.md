# `bayesTLS` R Package

<!-- badges: start -->
[![R-CMD-check](https://github.com/daniel1noble/bayesTLS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/daniel1noble/bayesTLS/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/daniel1noble/bayesTLS/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/daniel1noble/bayesTLS/actions/workflows/test-coverage.yaml)
[![Codecov test coverage](https://codecov.io/gh/daniel1noble/bayesTLS/graph/badge.svg?token=m2j2yA8jbt)](https://app.codecov.io/gh/daniel1noble/bayesTLS)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Ask Me Anything](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](https://github.com/daniel1noble/bayesTLS/issues)
[![Open Source](https://img.shields.io/badge/Open%20Source-%E2%9D%A4-brightgreen.svg)](https://github.com/daniel1noble/bayesTLS)
<!-- badges: end -->

`bayesTLS` is an R package for fitting joint Bayesian four-parameter logistic
(4PL) models to thermal death time (TDT) and thermal load sensitivity (TLS) data.
It works with raw counts or continuous proportions, extracts classical quantities
such as *z*, CTmax, *T*~crit~, and LT curves as posterior summaries, and provides
tools for uncertainty propagation, group contrasts, heat-injury prediction, and
bundled case-study data.

## Citation

Please cite the companion paper when using `bayesTLS`:

> Noble, D.W.A., Arnold, P.A., Nakagawa, S. & Pottier, P. (in preparation). A
> flexible modelling framework for estimating thermal tolerance and sensitivity.

This paper is the primary citation for both the statistical framework and the
`bayesTLS` package. Citation details will be updated here when the paper is
published. You can also run `citation("bayesTLS")` in R for the formatted
reference and a BibTeX entry.

## Installing the `bayesTLS` package

The analytical workflow ships as an installable R package at the root of this
repository. Install from GitHub:

```r
install.packages("remotes")  # if needed
remotes::install_github("daniel1noble/bayesTLS")
library(bayesTLS)
```

The package depends on [`brms`](https://paulbuerkner.com/brms/) (which needs a
Stan backend; [`cmdstanr`](https://mc-stan.org/cmdstanr/) is recommended). All
other imported dependencies (`dplyr`, `ggplot2`, `MASS`, `patchwork`,
`posterior`, `tibble`) are CRAN packages and resolve automatically. The classical
two-stage beta-binomial Stage 1 additionally uses the suggested
[`glmmTMB`](https://glmmtmb.github.io/glmmTMB/) package.

## Key functions in `bayesTLS`

The coding workflow is **decoupled** — each step is a standalone function so users can stop, swap, or skip at any stage:

| Function | Role |
|---|---|
| `standardize_data()` | Rename user columns to standard names; supports count responses (survived/dead/total) and continuous proportions; records response type + metadata. |
| `make_4pl_priors()` / `make_4pl_formula()` | Disjoint-bounds reparameterised 4PL; asymptote bounds adjustable for sublethal/PSII ranges. |
| `fit_4pl()` | Joint Bayesian 4PL via `brms`. Family resolved from the response — `beta_binomial(identity)` for counts, `Beta(identity)` for proportions. Centred temperature; `temp_effects` selects which 4PL parameters depend on temperature (all four by default; `"mid"` for sparse designs). |
| `extract_tdt()` | *z*, CTmax at the reference duration, and (lethal data only, `lethal = TRUE`) *T*~crit~ — all posterior transforms of the fit. *z* and CTmax share one posterior subsample. **Group-aware**: a fit whose CTmax/*z* vary by a moderator returns per-group results (auto-detected, or choose moderators with `by =`); a single-condition fit returns the usual one-group output. `seed =` makes the subsample + *T*~crit~ draws reproducible. |
| `tls()` / `tls_z()` / `tls_ctmax()` / `tls_tcrit()` | One call to derive *z*, CTmax, and (lethal) *T*~crit~ **per moderator group** (`by =`; sex, life stage, clone, …) from any fitted 4PL — a `bayes_tls` workflow *or* a hand-written `brms` 4PL — all from one consistent posterior. |
| `derive_z()`, `derive_tdt_curve()`, `derive_temperature_for_duration()`, `derive_tdt_landscape()` | Exported primitives behind `extract_tdt()`; `derive_z()` and `derive_temperature_for_duration()` are likewise group-aware (`by =`). |
| `extract_4pl_pars()`, `tdt_parameter_table()`, `tdt_quantile()`, `summarise_observed_survival()` | Per-draw 4PL parameters, natural-scale posterior parameter tables (group-aware via `by =`), TDT-friendly quantile summaries, and observed-survival cell summaries (mean ± SE per temp × duration). |
| `predict_survival_curves()`, `predict_heat_injury()`, `repair_rate_schoolfield()` | Posterior survival curves and heat-injury accumulation under a temperature trace (both per-group via `by =`), with an optional Sharpe–Schoolfield repair kernel supplied by `repair_rate_schoolfield()`. |
| `make_temperature_scenarios()`, `planted_dose_from_trace()` | Reference temperature traces and the analytical HI integral for validation. |
| `ts_stage1()`, `ts_stage2()`, `ts_ci()`, `ts_curve()` | The **classical two-stage TDT pipeline** for side-by-side comparison: per-temperature dose-response (`stats::glm` binomial or `glmmTMB` beta-binomial) → OLS of log10(LT50) on temperature → *z*, CTmax, *T*~crit~, with delta-method (Normal/*t*) or MVN intervals. One implementation, shared by the manuscript tutorial, case studies, and the bias simulation. |
| `diagnose_tdt_fit()` | Sampling diagnostics (R-hat, ESS, divergences) for a fitted workflow. |
| `plot_*()` + `theme_tdt()` | Plotting helpers (survival curves, TDT curve, tolerance landscape, heat injury, temperature scenarios/density, repair TPC) with a shared project theme. |
| `get_brmsfit()`, `get_z_draws()` / `get_z_summary()`, `get_ctmax_draws()` / `get_ctmax_summary()`, `get_tcrit_draws()` / `get_tcrit_summary()`, `get_surv_draws()`, `get_hi_draws()`, `has_fit()`, … | Accessors for the underlying fit, posterior draws, and posterior summaries. |
| `clock_to_minutes()`, `format_interval()` | Utilities: parse clock-time strings to elapsed minutes; format a posterior median with its interval as `median [lower, upper]`. |

Full reference: `?fit_4pl`, `?extract_tdt`, `?predict_heat_injury`, etc.

## `bayesTLS` Data

Five datasets are included with the `bayesTLS` package to make it easier to reproduce analyses and results in the paper and for testing purposes. These datasets are:

| Dataset | Help | Endpoint |
|---|---|---|
| `shrimp_lethal` | `?shrimp_lethal` | Brown-shrimp lethal TDT (survival counts) |
| `shrimp_sublethal` | `?shrimp_sublethal` | Brown-shrimp sublethal time-to-knockdown |
| `zebrafish_lethal` | `?zebrafish_lethal` | Zebrafish lethal TDT across life stages |
| `snowgum_psii` | `?snowgum_psii` | Snowgum leaf PSII (continuous proportion) |
| `dsuzukii` | `?dsuzukii` | *Drosophila suzukii* multi-trait TDT (lethal, knockdown, fertility; per individual) |

Each dataset can be loaded easily using `data(shrimp_lethal)` (as an example of loading the `shrimp_lethal` dataset). If you want other datasets loaded then simply replace `shrimp_lethal`. If you want to learn more about a dataset you can explore it's helpful `?shrimp_lethal`. 

## Reproducing the associated paper

This repository also contains the code needed to reproduce the companion paper,
its supplement, and the simulation results. To do this, install the package from
this repository, render the Quarto documents through the `Makefile`, and use the
tests to verify the package functions used by the manuscript and supplement. All
code chunks are provided in `ms/ms.qmd` and `ms/supplement.qmd`; use `make supp`
to render the supplement and `scripts/simulations/run_simulations.R` to rerun the
full two-stage-bias simulation suite. Run all scenarios with
`Rscript scripts/simulations/run_simulations.R`, or a subset by passing scenario
labels, e.g. `Rscript scripts/simulations/run_simulations.R scen9_tmax_060`. Each
scenario and the per-simulation pipeline are visible in that one script; the
reusable building blocks live in `scripts/simulations/sim_functions.R`.

### Repository layout

Use this map to find the files and outputs involved in reproduction:

- `Makefile` — main entry point for rendering the manuscript and supplement.
- `ms/ms.qmd` — manuscript source file; renders to `_output/ms.html`,
  `_output/ms.docx`, and `_output/ms.pdf`.
- `ms/supplement.qmd` — supplement source file; renders to `_output/supp.html`,
  `_output/supp.docx`, and `_output/supp.pdf`.
- `R/` — package source for the analysis functions used by the paper; see
  [Key functions](#key-functions-in-bayestls).
- `data/` — bundled package datasets as `.rda` files.
- `inst/extdata/` — source CSV files used to build the bundled datasets.
- `data-raw/make_datasets.R` — script that rebuilds the `.rda` datasets from the
  source CSV files.
- `output/models/` — cached `brms` fits (`.rds`) created and reused during
  renders.
- `scripts/simulations/` — the two-stage-bias simulation: `run_simulations.R`
  (scenarios + runner) and `sim_functions.R` (reusable building blocks).
- `scripts/` — other support scripts (presentation figures, document cleanup).
- `tests/testthat/` — fast unit tests plus optional `brms` integration tests.
- `bib/` — bibliography and the Ecology Letters CSL file.
- `notes/` — dated development notes for derivations, literature checks, and
  simulations.
- `man/` — roxygen-generated package help files.
- `_output/` — final rendered manuscript and supplement files; created by the
  render commands below.

### Reproducing the analysis

1. Install the package from this repository so the renders use the same code you
   are reproducing:

   ```r
   # install.packages("devtools")  # if needed
   devtools::install()
   ```

2. Render the manuscript, supplement, or both from the repository root:

   ```sh
   make all      # manuscript and supplement, all formats
   make ms       # manuscript only: _output/ms.{html,docx,pdf}
   make supp     # supplement only: _output/supp.{html,docx,pdf}
   ```

   You can also render one format at a time, for example `make ms-pdf` or
   `make supp-html`.

   The first supplement render may take several minutes because it fits or checks
   cached models. Subsequent renders reuse the `.rds` files in `output/models/`
   and should be much faster.

3. Inspect the rendered files in `_output/`. Use `make clean` to remove only
   `_output/`, or `make build-clean` to also remove the out-of-tree Quarto build
   cache at `~/Library/Caches/tls-render/`.

4. Run the package tests if you want to verify the analysis functions separately
   from the manuscript render:

   ```r
   devtools::test()

   Sys.setenv(RUN_BRMS_TESTS = "true")
   devtools::test()
   ```

   The first command runs the fast unit tests. The second also runs the gated
   `brms` integration tests, which fit small cached models and check parameter
   recovery.

### Rendering details

The manuscript and supplement are independent Quarto documents. There is no
project-level `_quarto.yml`; each `.qmd` file carries its own YAML so the files
can be rendered separately. The `Makefile` renders outside the Dropbox-synced
project directory, then copies only the final HTML, DOCX, and PDF files back to
`_output/`. This avoids conflicted intermediate files while preserving cached
model fits in `output/models/`.

Within a document, use Quarto's `@` references as normal (`@fig-X`, `@tbl-X`,
`@eq-X`, `@sec-X`). Between the manuscript and supplement, use plain text, such
as "Equation 7 of the manuscript" in the supplement or "Figure S2" in the
manuscript. The supplement render automatically labels figures, tables, and
equations with an "S" prefix.

See [CLAUDE.md §8a](CLAUDE.md) for more detail on the render architecture,
cross-reference handling, and the `brms` cache.

### Case studies in the supplement

The supplement works through the same workflow on simulated data and four
empirical examples:

1. **Brown shrimp** — lethal and sublethal endpoints.
2. **Zebrafish** — lethal TDT across life stages, fit as separate per-stage 4PLs
   and as a joint model with `life_stage` as a covariate.
3. **Snow gum leaf PSII** — continuous-proportion response fit with a Beta
   likelihood.
4. **Drosophila suzukii** — multi-trait TDT data, including mortality,
   knockdown, and fertility endpoints.

The supplement also contains the two-stage-bias simulation results and
sensitivity sweeps used to check the modelling framework.
