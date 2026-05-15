# bayesTLS

Methods paper on increasing the statistical rigor of the thermal death time (TDT) framework for ectotherm thermal tolerance. Target journal: *Methods in Ecology and Evolution*. The companion R package, [`bayesTLS`](https://github.com/daniel1noble/bayesTLS), ships from this repository.

## Purpose

The TDT framework is dominant in thermal biology but rests on a two-stage pipeline — per-temperature LT50 extraction followed by OLS regression on `log10(LT50)` against assay temperature — that violates several standard statistical assumptions: generated regressors, heteroscedastic point estimates, censored observations dropped or capped, proportion data handled with obsolete transformations, unmodelled overdispersion, and boundary observations treated by ad hoc fixes. This project catalogues those problems and demonstrates — on a zebrafish survival dataset — that a single joint Bayesian hierarchical model fit to the raw count data recovers the TDT / CTmax quantities of interest with properly propagated uncertainty.

The full problem catalogue with references lives in [tdt_problems.qmd](tdt_problems.qmd). Read it before contributing.

## Installing the `bayesTLS` R package

The analytical workflow is shipped as an installable R package at the root of this repository. Install from GitHub:

```r
# install.packages("remotes")  # if needed
remotes::install_github("daniel1noble/bayesTLS")
library(bayesTLS)
```

The package depends on [`brms`](https://paulbuerkner.com/brms/) (which in turn needs a Stan backend; [`cmdstanr`](https://mc-stan.org/cmdstanr/) is the recommended one). All other dependencies (`dplyr`, `ggplot2`, `patchwork`, `posterior`, `tibble`) are CRAN packages and resolve automatically.

Function reference is in the package documentation — e.g. `?fit_4pl`, `?extract_tdt`, `?predict_heat_injury`. The full guided walkthrough with simulated data plus the brown-shrimp case study lives in [`ms/supplement.qmd`](ms/supplement.qmd) and the rendered outputs in `_output/`.

## Reproducing the analysis

1. Install the package as above.
2. Render the supplement with `make supp` (HTML/DOCX/PDF) — the simulation tutorial fits and the shrimp case-study fit are cached to `output/models/` via `brms`'s `file_refit = "on_change"` mechanism, so first render takes a few minutes and subsequent renders are near-instant.
3. Inspect tests with `devtools::test()` (fast unit tests) or `RUN_BRMS_TESTS=true devtools::test()` (full integration suite that fits a small cached model and checks recovery).

## Rendering the manuscript and supplement

Two source documents, each rendered to HTML, DOCX, and PDF:

- `ms/ms.qmd` — manuscript
- `ms/supplement.qmd` — supplement

To produce all six outputs:

```sh
make all
```

Or one document at a time:

```sh
make ms     # manuscript only — _output/ms.{html,docx,pdf}
make supp   # supplement only — _output/supp.{html,docx,pdf}
```

Or one format at a time:

```sh
make ms-pdf      make ms-html      make ms-docx
make supp-pdf    make supp-html    make supp-docx
```

Final outputs land in `_output/`. The first render takes a few minutes (LaTeX, knitr cache); later renders reuse caches and are faster. Run `make clean` to remove `_output/`, or `make build-clean` to also wipe the out-of-tree build cache at `~/Library/Caches/tls-render/`.

**Cross-references.** Within a document, use Quarto's `@`-refs as normal (`@fig-X`, `@tbl-X`, `@eq-X`, `@sec-X`). Between the two documents, use plain text — e.g. write *"Equation 7 of the manuscript"* in the supplement, or *"Figure S2 in the supplement"* in the manuscript. The supplement labels figures/tables/equations with an "S" prefix automatically.

See [CLAUDE.md §8a](CLAUDE.md) for the full render architecture (out-of-tree builds, "S" label mechanism, brms cache via `output/models/`, Dropbox xattr setup).

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

## Planning log

### 2026-05-12 — Architecture & scope

**Summary.** Co-author planning discussion. Locked the scope of this paper, sketched a decoupled function pipeline, and split the static-CTmax-to-dynamic-CTmax conversion problem out into a future companion paper.

**Scope of this paper.**
- Joint Bayesian 4PL on proportion-style TDT data (survived / dead / total counts, durations, optional random effects).
- Sublethal / time-to-event data (knockdown time, fertility, photosystem-II failure) handled via the linear `log10(time) ~ temperature` route — z = -1/slope, CTmax_τ read off directly. These datasets often need heteroscedastic likelihoods because variance grows with temperature.
- Heat-injury accumulation under fluctuating temperature regimes, with and without repair.
- *Out of scope* for this paper: converting static CTmax measurements (single ramping rate, one CTmax) back to TDT parameters. That becomes Paper 3.

**Three-paper roadmap.**
1. *This paper.* TDT framework via joint Bayesian 4PL on proportion data + linear TDT for time-to-event/sublethal data + heat injury/repair.

**Different Paper**
2. Static-CTmax-to-TDT conversion. Likelihood- or simulation-based inference for the distribution of z that would have produced an observed CTmax × ramping-rate pattern, with uncertainty propagation. Empirical validation across datasets where both pipelines exist (zebrafish embryo, shrimp, our own). Demonstrates that the current literature's fixed-z assumption is unsupportable — across the 5 papers checked, z ranges ~1–7; sublethal-assay z is typically 1–2; assuming a single value silently changes conclusions.

**Pipeline architecture — decoupled, step-wise functions.** Each step is a standalone function so users can stop, swap, or skip:

1. `standardize_data()` — rename user columns to project-standard names (survived / dead / total / proportion, time, duration, optional random effects, optional covariates). Standardises fixed-effect naming so every downstream function reads the same object. Already partially in place; just needs polish and docs.
2. `make_4pl_priors()` — proportion data; bounds 0–1.
3. `make_TDT_priors()` — sublethal / unbounded data; bounds shifted to the data scale (e.g. photosystem II survives between 0.85 and 1).
4. `fit_4pl()` — joint Bayesian fit, beta-binomial likelihood. Centred temperature by default (decorrelates intercept/slope, fits where the data is, avoids extrapolation to T = 0). Default model = temperature on all four parameters (ℓ, u, k, mid); if a parameter has no temperature effect the term shrinks to zero and the unified correction factor degrades cleanly to the simpler case. Always center/ use Tref at 60 min by default but can be overwritten if needed. **Pending check: overfitting risk on small datasets — test before adopting.** On exit, print a console message summarising the assumed model, priors, bounds, and how to override.
5. `extract_tdt(fit, T_ref, target_surv = 0.5, TC_thresh = 0.05)` — z = -1/slope; CTmax (temperature at `target_surv` survival, default 0.5 → LT50); **T_crit** (temperature at `1 - TC_thresh` survival, default 0.05 → LD5). Same machinery for both: solve the 4PL inverse for temperature at the reference exposure time. Uses population-level estimates only; random effects irrelevant for extraction. Returns per-curve CTmax, per-curve T_crit, overall combined z, and full posterior draws for all three.
6. `predict_survival_curves()` — posterior predictive plots for sanity-checking the fit.
7. `predict_heat_injury(temp_series, z_post, CTmax_post, T_ref, repair = FALSE)` — accumulate damage under a temperature trace. Toggle returns heat-injury or survival. `repair = TRUE` uses Sharpe-Schoolfield with defaults estimated from where exponential damage accumulation begins (and an optimum temperature for repair sitting X°C below). Use Kelvin throughout the Arrhenius parts.
8. `fit_linear_TDT()` + matching `extract_tdt_linear()` for the time-to-event branch. Heteroscedastic likelihood option. Outputs match `extract_tdt()` so downstream functions don't care which path produced z and CTmax.
9. Plotting helpers for each of the above; faceting when interactions/categories are present. Plotting the TDT lines, time ~ temperature, on linear and log scale. Survival curves for plotting and then TDT landscape. 

**Supplement tutorial structure.** Intro → list of all exported functions with one-line descriptions → simulation walk-through that mirrors the paper section by section → final section showing the more complex case (K varies with T, additional covariates, interactions) using the same function set.

**Validation.** Heat-injury simulations: a flat trace (no injury), a single spike sized to deliver 50% mortality / 100% HI, a multi-spike trace. Cross-check `predict_heat_injury` recovers the planted dose under each.

**Centred temperature is the standard.** `T_ref` is the 60 min assay temperature used to centre. All downstream interpretation (CTmax at T_ref, z, heat injury) flows from this choice — document it prominently. 

### 2026-05-12 (cont'd) — `T_crit` lives in `extract_tdt()`

**Decision.** `T_crit` is added as a third quantity returned by `extract_tdt()`, alongside z and CTmax. New argument `TC_thresh` (default 0.05) controls the mortality threshold that defines it.

**Mechanically the same as CTmax.** CTmax is the temperature at which `target_surv = 0.5` (i.e. LT50) is reached at the reference exposure time. T_crit is the temperature at which `1 - TC_thresh = 0.95` survival is reached at the reference exposure time — i.e. the LD5 temperature. Both are the 4PL inverse solved for temperature at a fixed time; only the survival target differs. One function, one inversion routine, two outputs.

**Why not in `predict_heat_injury()`.** Putting T_crit there would force users to supply a temperature trace just to read off a quantity that is intrinsic to the dose-response surface. Many users will only have the fit. `extract_tdt()` is the right home — it operates on the fit alone.

**Why 5% (default).** 1% is within sampling randomness; 5% is robust and roughly matches where heat-injury accumulation visibly lifts off zero in `predict_heat_injury()` output. Threshold is exposed (`TC_thresh`) so users can tighten or loosen. Worth a visual sanity check against the heat-injury plots once the function is in place.

**Repair and T_crit.** Conceptually T_crit is where damage rate equals repair rate. But the TDT surface already encodes the *net* organismal response — repair is baked in. So T_crit can be defined from the TDT fit alone; the repair function (`predict_heat_injury(..., repair = TRUE)`) is for forward predictions under fluctuating temperatures, not for defining T_crit.

**Posterior T_crit — a contribution in its own right.** T_crit has historically been treated as a poorly-defined point value. Here it inherits the model's posterior — full uncertainty distribution, naturally wider at the low-mortality (small-effect, sparse-data) end. Worth emphasising in the manuscript: T_crit comes out *with* a credible interval, from the same single fit that gives z and CTmax. No back-and-forth, no separate experiment.

### 2026-05-12 (cont'd) — `T_crit` redefined as rate-multiplier (Faber et al. 2026)

**Decision.** Replaced the LD5-at-1-hour `T_crit` with the rate-multiplier definition from Faber et al. (2026, EcoEvoRxiv preprint `faber2026`, DOI: [10.32942/X2SM1B](https://doi.org/10.32942/X2SM1B)). For each posterior draw, sample $r^*$ uniformly on $\log_{10}$ over $[0.1, 1]$ % HI per hour, then $T_{crit} = CT_{max,1hr} + z \cdot \log_{10}(r^*/100)$. Pooled posterior carries both parameter uncertainty (in $CT_{max,1hr}$, $z$) and operational uncertainty (in $r^*$). Median sits at $CT_{max,1hr} - 2.5z$. Verified on the shrimp data — recovers the hand-picked $T_c = 25\,^\circ\text{C}$ used in the prior pipeline.

### 2026-05-13 — Case studies, conceptual figure, simulation, conclusions

**Case studies — four total, integrated rather than per-section.** Brown shrimp (lethal + sublethal — already in supplement); zebrafish (Pete is integrating now, using `fit_4pl` with `beta_binomial(link = "identity")`); a plant species with preprint precedence (z already known from prior work, so it acts as additional validation); a fourth plant species (TBD). Adrian/students have ~100 more species coming from a Southwest-China transplant project — useful for future validation but out of scope for this paper.

**Presentation format.** Single short manuscript section introducing the set with a one-paragraph blurb per case study (trait measured, experimental design, citation). Heavy lifting (data prep, fits, diagnostics, posteriors) lives in the supplement, one supplement subsection per case study. The manuscript carries the *integrated* figures only.

**Manuscript figures — three integrated figures, no per-case-study figures in ms.**

1. **Figure 1 — Conceptual figure (3 panels).** Panel A: per-temperature dose-response curves (two-stage Stage 1). Panel B: $\log_{10}\text{LT50}$ regression on temperature (two-stage Stage 2). Panel C: joint Bayesian 4PL — sketched two ways to be settled on a whiteboard:
   - *Option A*: same dose-response/LT50 view as A+B but as a single coupled model with credible bands on the LT50 line ($z$ = slope, $CT_{max,1hr}$ = where the line crosses $\log_{10} 60$).
   - *Option B*: the survival landscape ($T \times t$) with the 50% isocline; $CT_{max,1hr}$ sits at the intersection of that isocline with the 1-hour duration; uncertainty shown as dashed contour ridges. More information-dense but further from how field biologists think about TDT.
2. **Figure 2 — z and CTmax distributions across case studies.** Two side-by-side density panels (z left, $CT_{max,1hr}$ right). For each case study, posterior density with 95% CrI bar; two-stage point estimate overlaid as a vertical line for direct comparison. Species silhouettes/iconography to make it visually appealing. Communicates the headline equivalence claim: joint 4PL and two-stage land in the same place on average, but the joint model gives a calibrated CrI for every case study.
3. **Figure 3 — Heat-injury and survival under one simulated temperature trace, applied to all 4 case studies.** A single simulated time series at the top (a realistic regime with one or two heat-wave events), then a row of $\sim 4$ panels — one per case study — showing posterior HI accumulation and predicted survival. Demonstrates the ecological interpretation: thermally sensitive species accumulate damage rapidly and lose population fraction, less sensitive species stay stable. Concrete forward-prediction utility. Conceptually mirrors the tree-paper figure with heat-wave spikes and species-specific responses.

**Two-stage bias simulation — new R script.** Run an explicit simulation harness that varies (a) sample size per cell, (b) overdispersion, (c) whether the assay design covers the full $[0, 1]$ survival range or only part of it, then compare joint Bayesian 4PL posterior estimates to two-stage point estimates of $z$ and $CT_{max,1hr}$ against the simulation truth. Goal: quantify the direction and magnitude of two-stage bias as a function of sample size and design coverage. Daniel's prior: bias is real and worst under small-sample / partial-coverage designs, but the magnitude is unknown. Outputs: bias plots, RMSE comparison, coverage of nominal 95% intervals. Lives in `scripts/sim_twostage_bias.R`; summary tables and one figure folded into the supplement.

**Conclusions / discussion blocks.**
- Joint 4PL recovers what two-stage already recovers (on average) — **enhancement, not replacement**.
- Sampling error becomes available for every quantity → makes meta-analyses tractable (point estimates *with* CrIs from every contributing study, in a common format).
- Interpretation is more user-friendly: predicted survival probabilities and HI trajectories speak to ecological audiences, where $z$ and $CT_{max,1hr}$ on their own do not.
- Uncertainty propagates to nature: forward predictions under measured/projected field temperature traces inherit the full posterior, which previous TDT pipelines never did.
- Future enhancements (out of scope for this paper, but worth pointing to): acclimation responses (Arnold; chat with Johannes in France), life-stage-specific fits, repair-function calibration.
- *Discount-data advantage*: two-stage TDT requires assay designs that fully cross the $[0, 1]$ survival range — partial-response data, control temperatures, and replicates that didn't fully kill the organism are typically discarded. The joint 4PL uses all of it; even constant-survival data at control temperatures help estimate the upper asymptote $u$.
- *Faber-et-al T_crit*: cross-method agreement (~1 °C across performance-assay, constant-TDT breakpoint, and alternating-CT methods in two species) suggests `T_crit` is reasonably well-identified empirically — but the uncertainty *between* methods sits on top of within-method sampling uncertainty, which we propagate explicitly.
- *Limitations*: assumptions about how the HI integral and repair are formulated; rate-multiplier $T_{crit}$ uses an operationally chosen $r^*$ range; the framework still relies on the 4PL functional form being adequate for the dose-response.

**Manuscript figure-1 conceptual sketch — open questions to decide on whiteboard.**
- How to draw the joint-model panel so it communicates *one fit, two quantities* (z + CTmax) rather than feeling like a third side-by-side procedure.
- Whether to fold $T_{crit}$ into the same figure or leave it as a downstream consequence with its own visualisation in the body.
- Whether the uncertainty should be drawn as ribbons on the LT50 line or as dashed contours on the landscape.

**Workflow notes (logistics).**
- Each new case study is best developed in a self-contained script (`scripts/case_study_<species>.R`) and then merged into the supplement once the fit is sane. This keeps the supplement clean while Pete iterates on the zebrafish and plant fits.
- The supplement is the HTML primary artefact (updatable post-publication via Zenodo/GitHub release); a PDF copy is fine but secondary.
- Repo is private. New collaborators need to push at least one commit (e.g., touch a file, `git push`) before `remotes::install_github()` resolves their identity.

## TODO list

### Completed (locked-in as of 2026-05-13)

**`bayesTLS` R package — installable from GitHub, R CMD check passing:**

- [x] `standardize_data()` — column standardisation, random-effect grouping, mean-centred temperature, attached metadata.
- [x] `make_4pl_priors()` — disjoint-bounds priors for proportion data; bounds adjustable via `lower` / `upper` for sublethal/PSII-style data.
- [x] `fit_4pl()` — joint Bayesian wrapper with `beta_binomial(link = "identity")` (default; `binomial(link = "identity")` available for no-overdispersion), disjoint-bounds reparameterisation, default `~ temp_c` slope on all four 4PL sub-parameters, mean-centred temperature, random intercepts on `mid` via `random_effects = c(...)`, brms `file_refit = "on_change"` caching.
- [x] `extract_tdt()` — z (per-draw OLS of $\log_{10}\text{LT50}$ on T), CTmax (per-draw 4PL inversion at $t_\text{ref}$), T_crit (rate-multiplier definition from Faber et al.: $T_{crit} = CT_{max,1hr} + z \cdot \log_{10}(r^*/100)$, $r^*$ sampled uniformly on $\log_{10}$ across `TC_rate_range`, default $[0.1, 1]$ %/hr).
- [x] `derive_tdt_curve()`, `derive_temperature_for_duration()`, `derive_tdt_parameters()`, `derive_tdt_landscape()` — the primitives behind `extract_tdt()`, exported so users can call them directly.
- [x] `predict_survival_curves()` — posterior survival curves on a (temp × duration) grid.
- [x] `predict_heat_injury()` — HI accumulation under any temperature trace; optional Sharpe-Schoolfield repair (`repair_rate_schoolfield()`), Kelvin internally, irreversible-mortality option, save-draws option.
- [x] `make_temperature_scenarios()` — four reference traces (flat / single-spike / multi-spike / diurnal). Single-spike is calibrated so a 1-hour spike at `CTmax_1hr` delivers ~100% LT50 dose by construction. Diurnal is multi-day with day-to-day variability in peak temperature.
- [x] `planted_dose_from_trace()` — analytical HI integral truth for validation.
- [x] Plotting helpers — `plot_survival_curves()` (default linear time + viridis temperature; classical log-time via `log_time = TRUE`), `plot_tdt_curve()`, `plot_tdt_landscape()` (default linear time), `plot_temperature_density()`, `plot_temperature_scenarios()`, `plot_heat_injury()`, `plot_repair_rate()`. Shared `theme_tdt()` for visual consistency.

**Validation & supplement content:**

- [x] Heat-injury simulation harness — three traces (flat / single spike / multi-spike), planted-dose validation. Single-spike at posterior median $CT_{max,1hr}$ recovers ~100% HI by construction (used as a self-check in the shrimp section).
- [x] Overfitting check for default `fit_4pl()` — supplement § *A worked example: temperature effects on every shape parameter* shows the fit picks up real T effects on $u$ and $k$ and shrinks them toward zero when absent.
- [x] Supplement tutorial — intro, function tour, simulation walk-through, joint Bayesian 4PL fit, $z$ / $CT_{max,1hr}$ / $T_{crit}$ derivation, comparison to the classical two-stage pipeline, T-varying-shape worked example, heat injury under reference + diurnal traces (with and without repair).
- [x] Sublethal / linear-fit pipeline — brown-shrimp sublethal time-to-knockdown analysed both via hierarchical linear Bayesian model (`log10_time ~ temp_c + (1|date) + (1|tank) + (1|cup)`) and via Ørsted-style 4PL on proportion-knocked-down data; side-by-side comparison plot and table in the supplement.

### This paper — outstanding

**Case studies — four total, integrated rather than per-section in the manuscript:**

- [x] Brown shrimp (lethal + sublethal) — already in the supplement.
- [ ] Zebrafish lethal TDT (Pete is integrating). Fit with `fit_4pl()` + `beta_binomial(link = "identity")`; re-use the plotting helpers. Develop in `scripts/case_study_zebrafish.R`, then merge into the supplement once the fit is sane.
- [ ] Plant case study #1 — Pete's species with preprint precedence (z already published; acts as additional cross-method validation).
- [ ] Plant case study #2 — Pete's second plant species (TBD; one clean dataset, not the broader ~100-species China-transplant set).
- [ ] Each case study contributes one supplement subsection (data prep → fit → diagnostics → `extract_tdt()` summary) and one row's worth of data to the integrated manuscript figures below.

**Manuscript figures — three integrated figures, no per-case-study figures in the manuscript:**

- [ ] **Figure 1 — Conceptual figure (3 panels).** Panel A: per-temperature dose-response curves (two-stage Stage 1). Panel B: $\log_{10}\text{LT50}$ regression on temperature (two-stage Stage 2). Panel C: joint Bayesian 4PL — settle the visual on whiteboard. Candidate views: (a) LT50-on-temperature line with CrI ribbons, matching A+B; (b) survival landscape with 50%-survival isocline and dashed CrI contour ridges. Decide which.
- [ ] **Figure 2 — Distribution of $z$ and $CT_{max,1hr}$ across case studies.** Two side-by-side density panels (z left, $CT_{max,1hr}$ right). One density per case study (4 species), 95% CrI bars, two-stage point estimates overlaid as vertical lines for direct comparison. Species silhouettes for visual appeal. Headline: joint 4PL and two-stage agree on average, but only the joint 4PL gives a calibrated CrI on every quantity.
- [ ] **Figure 3 — Heat injury + survival under one simulated temperature trace, applied to all 4 case studies.** Top: a single simulated time series (a realistic regime with one or two heat-wave events, calibrated so some species accumulate damage and others stay stable). Below: 4 panels (one per species) showing posterior HI accumulation and predicted survival. Demonstrates the ecological/community-level interpretation of the framework.

**Two-stage bias simulation:**

- [ ] `scripts/sim_twostage_bias.R` — self-contained harness varying (a) sample size per cell (small / moderate / large), (b) overdispersion ($\phi$), (c) design coverage (full $[0, 1]$ survival range vs partial). For each simulated dataset, fit both the joint Bayesian 4PL and the two-stage pipeline; report bias, RMSE, and 95% interval coverage for $z$ and $CT_{max,1hr}$ against the simulation truth.
- [ ] Save outputs (per-condition summaries) to `output/sim_twostage/` so the supplement can load and summarise them without re-running the harness at render time.
- [ ] Fold one figure and one or two summary tables into the supplement.

**Manuscript prose (after the three figures are sketched):**

- [ ] Single case-study intro section in `ms/ms.qmd` — brief paragraph per species (trait measured, experimental design, n temperatures, n durations, citation). Heavy lifting stays in the supplement.
- [ ] Discussion / conclusions — pull from the 2026-05-13 planning log entry above:
  - Equivalence-as-enhancement framing (joint 4PL reproduces what two-stage already gives, on average).
  - Sampling error enables meta-analysis (point estimates *with* CrIs, in a common format).
  - Interpretation gains (survival probabilities + HI trajectories speak to ecological audiences).
  - Uncertainty propagation to forward predictions in nature.
  - Discount-data advantage (joint 4PL uses partial-response and control data the two-stage discards).
  - Faber-et-al cross-method $T_{crit}$ agreement vs propagated within-method uncertainty.
- [ ] Limitations — HI integral assumptions, $r^*$ range choice for $T_{crit}$, 4PL functional-form adequacy. Be explicit but not apologetic.
- [ ] Future directions (signposted, not implemented in this paper): acclimation responses (Arnold; Johannes chat in France), life-stage-specific fits, repair-function calibration from sub-lethal panels rather than fixed Sharpe-Schoolfield parameters.

**Reference housekeeping:**

- [ ] Pete to dump any of his references that fit naturally into the existing prose.
- [ ] Daniel to chase down the references flagged as "book — needs primary-source verification" and replace with primary citations where possible.

### Deferred (`fit_linear_TDT`)

The sublethal section of the supplement currently calls `brms::brm()` inline for the linear time-to-event model. The original architecture planned a `fit_linear_TDT()` wrapper exposing the same downstream interface (z, CTmax, posterior draws) as `fit_4pl()`. Useful for API consistency; not blocking this paper. Schedule for a post-submission package release.

- [ ] `fit_linear_TDT()` + `extract_tdt_linear()` — same downstream interface as `fit_4pl()` / `extract_tdt()`; heteroscedastic likelihood option.

### Deferred to companion papers

- Paper 2 — flow-cytometry / cell-level application of the framework.
- Paper 3 — static-CTmax-to-TDT conversion.
  - [ ] (later) Likelihood/simulation function to infer z-distribution from CTmax × ramping-rate data with uncertainty propagation.
  - [ ] (later) Validation across datasets with both pipelines.
  - [ ] (later) Show fixed-z assumption (e.g. Jørgensen 2021) is unsupportable given observed z variation.





