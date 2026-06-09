# bayesTLS

<!-- badges: start -->
[![R-CMD-check](https://github.com/daniel1noble/bayesTLS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/daniel1noble/bayesTLS/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/daniel1noble/bayesTLS/graph/badge.svg)](https://app.codecov.io/gh/daniel1noble/bayesTLS)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
<!-- badges: end -->

A methods project increasing the statistical rigour of the thermal death time
(TDT) / thermal load sensitivity framework for ectotherm thermal tolerance. The
repository is simultaneously (a) the installable R package
[`bayesTLS`](https://github.com/daniel1noble/bayesTLS) and (b) the manuscript
compendium (`ms/`, `notes/`, bundled data). Target journal: *Ecology Letters*
(Methods).

## Purpose

The classical TDT pipeline is two-stage: per-temperature LT50 point estimates are
extracted first, then OLS is fit to `log10(LT50)` against assay temperature. This
loses uncertainty (Stage-1 estimates enter Stage-2 as if known), mishandles
proportion/count data, drops or caps censored observations, ignores
overdispersion and hierarchical structure, and hides analyst degrees of freedom.

We replace it with a single joint Bayesian hierarchical four-parameter logistic
(4PL) model fit to the raw counts (or continuous proportions). The classical
quantities — *z*, CTmax at a reference duration, *T*~crit~, LT curves,
heat-injury trajectories — are recovered as posterior summaries of that one fit,
with uncertainty propagated throughout. The framing is *equivalence-as-enhancement*:
under benign conditions the joint 4PL reproduces the field's existing point
estimates, while additionally delivering calibrated intervals on every derived
quantity, group contrasts, and forward heat-injury prediction that the two-stage
pipeline cannot.

The full problem catalogue with references is in
[notes/tdt_problems.qmd](notes/tdt_problems.qmd). Read it before contributing.

## Repository layout

The repo follows standard R-package layout *and* carries the manuscript
compendium:

- `R/` — package source (the analysis pipeline; see [Key functions](#key-functions)).
- `tests/testthat/` — unit tests plus a gated brms integration suite.
- `data/` — bundled datasets as `.rda` (the five case-study datasets); built from
  `inst/extdata/` CSVs by `data-raw/make_datasets.R`.
- `inst/extdata/` — raw CSVs shipped with the package.
- `man/` — roxygen-generated help.
- `ms/ms.qmd`, `ms/supplement.qmd` — manuscript and supplement, each
  self-contained (there is **no** project-level `_quarto.yml`).
- `notes/` — dated `.qmd` notes, one per compartmentalised problem (derivations,
  literature checks, simulations). `notes/legacy/` holds superseded prototype
  notebooks.
- `output/models/` — cached `brms` fits (`.rds`), reused across renders.
- `bib/` — `tdt_problems.bib` (auto-exported from Zotero — never edit directly)
  plus the Ecology Letters CSL.
- `scripts/` — the two-stage-bias simulation harness and data-prep helpers.

## Installing the `bayesTLS` package

The analytical workflow ships as an installable R package at the root of this
repository. Install from GitHub:

```r
# install.packages("remotes")  # if needed
remotes::install_github("daniel1noble/bayesTLS")
library(bayesTLS)
```

The package depends on [`brms`](https://paulbuerkner.com/brms/) (which needs a
Stan backend; [`cmdstanr`](https://mc-stan.org/cmdstanr/) is recommended). All
other dependencies (`dplyr`, `ggplot2`, `patchwork`, `posterior`, `tibble`) are
CRAN packages and resolve automatically.

## Key functions

The pipeline is **decoupled** — each step is a standalone function so users can
stop, swap, or skip at any stage:

| Function | Role |
|---|---|
| `standardize_data()` | Rename user columns to standard names; supports count responses (survived/dead/total) and continuous proportions; records response type + metadata. |
| `make_4pl_priors()` / `make_4pl_formula()` | Disjoint-bounds reparameterised 4PL; asymptote bounds adjustable for sublethal/PSII ranges. |
| `fit_4pl()` | Joint Bayesian 4PL via `brms`. Family resolved from the response — `beta_binomial(identity)` for counts, `Beta(identity)` for proportions. Centred temperature; `temp_effects` selects which 4PL parameters depend on temperature (all four by default; `"mid"` for sparse designs). |
| `extract_tdt()` | *z*, CTmax at the reference duration, and (lethal data only, `lethal = TRUE`) *T*~crit~ — all posterior transforms of the single fit. *z* and CTmax share one posterior subsample. |
| `derive_z()`, `derive_tdt_curve()`, `derive_temperature_for_duration()`, `derive_tdt_landscape()` | Exported primitives behind `extract_tdt()`. |
| `predict_survival_curves()`, `predict_heat_injury()` | Posterior survival curves; heat-injury accumulation under a temperature trace (optional Sharpe-Schoolfield repair). |
| `make_temperature_scenarios()`, `planted_dose_from_trace()` | Reference temperature traces and the analytical HI integral for validation. |
| `plot_*()` + `theme_tdt()` | Plotting helpers with a shared project theme. |
| `get_brmsfit()`, `get_z_draws()`, `get_ctmax_draws()`, `get_tcrit_draws()`, … | Accessors for the underlying fit and posterior draws. |

Full reference: `?fit_4pl`, `?extract_tdt`, `?predict_heat_injury`, etc.

## Reproducing the analysis

1. Install the package as above.
2. Render the supplement with `make supp` — the simulation tutorial fits and the
   case-study fits are cached to `output/models/`, so the first render takes a
   few minutes and subsequent renders are near-instant.
3. Run the tests with `devtools::test()` (fast unit tests) or
   `RUN_BRMS_TESTS=true devtools::test()` (full integration suite that fits small
   cached models and checks parameter recovery).

## Rendering the manuscript and supplement

Two source documents, each rendered to HTML, DOCX, and PDF:

- `ms/ms.qmd` — manuscript
- `ms/supplement.qmd` — supplement

```sh
make all      # all six outputs
make ms       # manuscript only — _output/ms.{html,docx,pdf}
make supp     # supplement only — _output/supp.{html,docx,pdf}
```

Or one format at a time (`make ms-pdf`, `make supp-html`, …). Final outputs land
in `_output/`. `make clean` removes `_output/`; `make build-clean` also wipes the
out-of-tree build cache at `~/Library/Caches/tls-render/`.

**Cross-references.** Within a document, use Quarto's `@`-refs as normal
(`@fig-X`, `@tbl-X`, `@eq-X`, `@sec-X`). Between the two documents, use plain text
— e.g. *"Equation 7 of the manuscript"* in the supplement, or *"Figure S2"* in
the manuscript. The supplement labels figures/tables/equations with an "S" prefix
automatically.

See [CLAUDE.md §8a](CLAUDE.md) for the full render architecture (out-of-tree
builds to avoid Dropbox conflicts, the "S" label mechanism, brms cache via
`output/models/`).

## Data

Five datasets ship with the package:

| Dataset | Help | Endpoint |
|---|---|---|
| `shrimp_lethal` | `?shrimp_lethal` | Brown-shrimp lethal TDT (survival counts) |
| `shrimp_sublethal` | `?shrimp_sublethal` | Brown-shrimp sublethal time-to-knockdown |
| `zebrafish_lethal` | `?zebrafish_lethal` | Zebrafish lethal TDT across life stages |
| `snowgum_psii` | `?snowgum_psii` | Snowgum leaf PSII (continuous proportion) |

Each is built from a raw CSV in `inst/extdata/` by `data-raw/make_datasets.R` and
carries column-by-column roxygen documentation in `man/`. Raw data are read-only.

## Case studies (worked in the supplement)

1. **Brown shrimp** — lethal + sublethal endpoints.
2. **Zebrafish** — lethal TDT across life stages, fit two ways (separate per-stage
   4PLs vs a joint model with `life_stage` as a covariate), with a relative-vs-absolute
   survival-threshold sensitivity check.
3. **Snowgum leaf PSII** — continuous-proportion (Beta) response.

The supplement also contains the full two-stage-bias simulation results and
sensitivity sweeps.

## Associated publications

None yet. Target journal: *Ecology Letters* (Methods). Noble, Arnold & Pottier
(in preparation).

## Conventions

Coding, writing, testing, reproducibility, and collaboration conventions live in
[CLAUDE.md](CLAUDE.md).
</content>
</invoke>
