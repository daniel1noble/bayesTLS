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

## TODO list

**TODO list:** 

- [ ]  Make a general function to fit 4 PL model from survival data
    - [ ]  Fix issues with priors being too strong (i.e., they are currently derived from data)
    - [ ]  Make sure to use beta binomial implementation instead of OLRE
    - [ ]  Add options for covariates + interactions + random effects?? Or just leave the code so people can adapt easily?
    - [ ]  Add option for number of chains, iterations etc
    - [ ]  Spit summary (z, R2 etc) + potential problems

- [ ]  Make a general function to generate survival curves and check everything is ok. 

- [ ]  Make a general function to derive TDT curve from 4PL model
    - [ ]  Make an adjustment with different temperature targets (0.5 by default, but could be less)
     - [ ]  If there are multiple categories (interaction), then facet_wrap
   
- [ ]  Make a general function to derive TDT landscape 

- [ ]  Make a general function to derive temperatures tolerated for 1 hour from 4PL model

- [ ]  Make a general function to predict accumulation of heat injury in fluctuating temperature regimes
    - [ ]  Simulate temperatures with known predicted accumulation of injury. One with none; one with a single spike (e.g., predicted to result in 50% mortality or 100% injury); one with multiple spikes.
    - [ ]  Have an option for repair = TRUE
    - [ ]  By default, we can try to estimate the TPC for repair based on when damage starts to accumulate exponentially; and X degrees below this for the optimum temperature for repair? 
    - [ ]  Make sure everything is in Kelvin if we are using Arrhenius equations. 

- [ ]  Make a general function to plot predicted injury/repair at different temperatures (with uncertainty)

- [ ]  Make a general function to predict survival/mortality rate over the time series. 

- [ ]  Make a general function to derive TDT parameters from dynamic data; and dynamic data to TDT parameters.





