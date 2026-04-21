# tls_model_equivalence

Methods paper on increasing the statistical rigor of the thermal death time (TDT) framework for ectotherm thermal tolerance. Target journal: *Methods in Ecology and Evolution*.

## Purpose

The TDT framework is dominant in thermal biology but rests on a two-stage pipeline — per-temperature LT50 extraction followed by OLS regression on `log10(LT50)` against assay temperature — that violates several standard statistical assumptions: generated regressors, heteroscedastic point estimates, censored observations dropped or capped, proportion data handled with obsolete transformations, unmodelled overdispersion, and boundary observations treated by ad hoc fixes. This project catalogues those problems and demonstrates — on a zebrafish survival dataset — that a single joint Bayesian hierarchical model fit to the raw count data recovers the TDT / CTmax quantities of interest with properly propagated uncertainty.

The full problem catalogue with references lives in [tdt_problems.qmd](tdt_problems.qmd). Read it before contributing.

## Status

Early stage. Currently in the repo:

- [tdt_problems.qmd](tdt_problems.qmd) — statistical critique of the TDT literature, with references.
- [R/bayesian_4pl.R](R/bayesian_4pl.R) — prototype joint Bayesian 4-parameter logistic fit (`brms` / `cmdstanr`) on zebrafish survival counts: observation-level overdispersion, experiment-day random effects, life-stage-varying asymptotes, slope, and inflection. Demonstrates the proposed alternative to the two-stage TDT pipeline.
- [bib/](bib/) — `.bib` file for `tdt_problems.qmd`, the Ecology Letters CSL style, and [bib/templates/](bib/templates/) holding the docx template used to render `.qmd` files to `.docx` for collaborators.

Planned per [CLAUDE.md](CLAUDE.md) but not yet present: `data/`, documented and tested functions in `R/`, `tests/testthat/`, `ms/` (manuscript + SI), `notes/` (per-problem derivations), `output/`, `.github/workflows/`, `renv.lock`.

## Reproducing the analysis

To be written once the compendium scaffold is in place. The current prototype in [R/bayesian_4pl.R](R/bayesian_4pl.R) expects `data/data_zebrafish_TDT.xlsx` (sheet `LETHAL_TDT`), which is not yet in the repo, and it references plotting helpers (`base_theme`, `life_stage_cols`, …) that have not yet been extracted into `R/`.

## Data

See `data/README.md` (to be added) for column-by-column documentation of each raw data file, along with source, license, ethics, and sharing restrictions.

## Key functions

To be written once functions are extracted from [R/bayesian_4pl.R](R/bayesian_4pl.R) into documented, tested units under `R/`.

## Associated publications

None yet. Target journal: *Methods in Ecology and Evolution*.

## Supplementary material

None yet. Any GitHub Pages tutorials or supplementary case studies will be linked here as they are written.

## Conventions

Coding, writing, testing, and collaboration conventions live in [CLAUDE.md](CLAUDE.md).
