# Codex Project Memory

Last reviewed: 2026-05-25

Working memory for AI-agent sessions in this repository: what the project is, the
current architecture, the code that matters, and the rough edges to remember
before editing. `CLAUDE.md` is the authoritative source for conventions and
non-negotiables; this file is orientation, not rules.

## Project identity

- Repository: `bayesTLS` (https://github.com/daniel1noble/bayesTLS); the local
  Dropbox folder keeps the historical name `tls_model_equivalence/`.
- The repo is **both** an installable R package (`R/`, `man/`, `DESCRIPTION`,
  `NAMESPACE`, `tests/`, `data/`, `inst/extdata/`) **and** the manuscript
  compendium (`ms/`, `notes/`, `bib/`, `output/`, `scripts/`).
- Target journal: **Ecology Letters** (Methods article). (Earlier drafts aimed at
  Methods in Ecology and Evolution — that target is retired.)
- Authors: Daniel W. A. Noble, Pieter A. Arnold, Patrice Pottier.
- Empirical scope: three worked case studies — brown shrimp (lethal + sublethal),
  zebrafish (lethal across life stages), and snowgum leaf PSII (continuous
  proportion).

### Scientific pitch

The classical thermal death time / thermal load sensitivity pipeline estimates
per-temperature LT50 point estimates and then fits OLS to `log10(LT50) ~ T`. This
loses uncertainty, mishandles proportion/count data, drops or caps censored
temperatures, ignores overdispersion and hierarchy, and hides analyst degrees of
freedom. The replacement is a single joint Bayesian hierarchical 4PL model fit to
the raw survival counts (or continuous proportions). Classical outputs — *z*,
CTmax at a reference duration, *T*~crit~, LT curves, heat-injury trajectories —
are derived from posterior draws of that one fit, not estimated in a second
stage.

**Framing (important):** lead with *equivalence-as-enhancement*. Under benign
conditions the joint 4PL reproduces the field's existing point estimates
(validating the existing literature); the recommendation is nonetheless the joint
4PL, because it adds calibrated uncertainty on derived quantities, group
contrasts, hierarchical structure, censoring, and dynamic heat-injury prediction
the two-stage cannot. Do **not** soften this to "use either" or "if results agree,
use the two-stage". Do **not** strawman the two-stage by only simulating
conditions that break it.

## Engagement rules (from CLAUDE.md)

- Act as a biostatistician + biologist, co-author plus critical reviewer; disagree
  constructively and separate standard practice / recommendation / unsettled.
- No fabricated citations, DOIs, p-values, effect sizes, or methods.
- Every non-obvious decision is logged in a dated `notes/*.qmd`. Render notes after
  editing them and inspect the printed output.
- Raw data are read-only. Derived objects go to `output/`.
- Reusable logic lives in `R/` with roxygen docs + a runnable example, and a
  matching test in `tests/testthat/`.
- Manuscript numbers, tables, and figures are produced by code — none typed by
  hand.
- Report simulation/data results as flat numbers + confounders not yet ruled out;
  no advocacy adjectives.

## Architecture and build (current)

There is **no project-level `_quarto.yml`** and **no Quarto book**. The earlier
book setup (`_quarto.yml`, `_book/`, `index.qmd`, `suppfig-*`/`supptbl-*` custom
floats) is gone — do not reintroduce it. Each `.qmd` carries its complete YAML so
renders are fully independent.

Two deliverables, each in HTML/DOCX/PDF, via the `Makefile`:

- `make ms`   → `_output/ms.{html,docx,pdf}` (from `ms/ms.qmd`)
- `make supp` → `_output/supp.{html,docx,pdf}` (from `ms/supplement.qmd`)
- `make all`  → all six.

Key build mechanics (see `Makefile` + CLAUDE.md §8a):

- **Out-of-tree builds.** The project lives in Dropbox; renders run in
  `~/Library/Caches/tls-render/` to avoid "conflicted copy" duplicates. Sources
  are rsync'd in; `bib/`, `output/`, `data/` are symlinked; only final artefacts
  are copied back to `_output/`. `R/` is **not** symlinked — the supplement loads
  the analytical functions via `library(bayesTLS)`. Install the package once
  before rendering.
- **"S" labels.** Supplement figures/tables/equations get "Figure S1, Table S1, …"
  via `--metadata-file ms/_supp-overrides.yml` (HTML/DOCX) and a Lua filter via
  `ms/_supp-pdf-overrides.yml` (PDF), applied at the command line — kept out of
  `supplement.qmd`'s own YAML. A post-render step strips the U+00A0 nbsp Pandoc
  inserts between "Figure S" and the number.
- **Cross-document refs** use plain text ("Equation 7 of the manuscript"); `@`-refs
  are within-document only.
- **brms cache.** Fits are cached as `.rds` in `output/models/` and reused via
  `file_refit`. Simulation-tutorial fits use `file_refit = "never"` (delete the
  `.rds` to force a refit); case-study fits use the default `"on_change"`.
- **bib.** `bib/tdt_problems.bib` is **auto-exported from Zotero** (Better BibTeX
  "keep updated"). Never edit it directly — add citations on the Zotero side.

CI from CLAUDE.md §6 (testthat + render + lint on push/PR) is **not yet wired up**
— there is no `.github/workflows/` directory. Don't claim CI runs.

## The package (`R/`)

Decoupled, stepwise pipeline (see CLAUDE.md memory + `notes/`):

- `standardize_data()` — rename user columns to standard names; supports **count**
  (survived/dead/total) and **continuous-proportion** responses. Records
  `response_type` ("count" | "proportion"), `response_var`, centred temperature,
  bounds, and random-effect grouping in the workflow metadata. `n_total` is
  optional (proportion path).
- `make_4pl_priors()` / `make_4pl_formula()` — disjoint-bounds reparameterised
  4PL. `compute_4pl_bounds(lower, upper, pad = 0.001, gap = 0.002)` splits the
  response range at its midpoint into two disjoint `inv_logit` intervals (for
  proportion data: low ∈ (0.001, 0.499), up ∈ (0.501, 0.999)); the gap kills
  label-switching, the pad keeps asymptotes off the exact boundaries. Bounds shift
  for sublethal/PSII ranges (e.g. lower = 0.85).
- `fit_4pl()` — joint Bayesian 4PL via native `brms`. Default family resolved from
  `response_type`: `beta_binomial(link = "identity")` for counts,
  `Beta(link = "identity")` for proportions. **Identity link throughout**
  (validated indistinguishable from logit on both response types). Centred
  temperature. `temp_effects` (subset of `c("low","up","k","mid")`, default all
  four; `mid` always required) controls which 4PL sub-parameters depend on
  temperature — use `"mid"` (classical constant-shape TDT) for sparse designs to
  avoid the divergences/unidentified-*z* the all-four default produces there.
- `extract_tdt(wf, t_ref = 60, lethal = FALSE, target_surv = 0.5, …)` — returns
  *z* and CTmax at `t_ref` always; *T*~crit~ only when `lethal = TRUE`. It is a
  **posterior transform, not a two-stage regression**: it pulls the
  population-level 4PL coefficient draws once, subsamples once, and computes both
  *z* and CTmax from that single set so they share draws by construction.
  - *z*: relative threshold (default) → `-1 / b_mid_temp_c` per draw (exact);
    absolute threshold → per-draw local `z(T)` via finite difference of the
    closed-form LT curve.
  - CTmax: closed-form crossing of `log10 LT(T) = log10(t_ref)`.
  - *T*~crit~ (lethal only): rate-multiplier definition of Faber et al. (2026) —
    per draw, sample `r*` uniformly on log10 across `TC_rate_range` (default
    `c(0.1, 1)` % HI/hr) and set `T_crit = CTmax_1hr + z * log10(r*/100)`. Median
    ≈ `CTmax_1hr − 2.5·z`; CrI noticeably wider than CTmax — report the full
    interval. For **sublethal** endpoints the fitted *z* is a performance-reduction
    slope, so the rate-multiplier formula does not apply (`lethal = FALSE`).
- `derive_z()`, `derive_tdt_curve()`, `derive_temperature_for_duration()`,
  `derive_tdt_landscape()` — exported primitives.
- `predict_survival_curves()`, `predict_heat_injury()` (Kelvin internally;
  optional Sharpe-Schoolfield repair via `repair_rate_schoolfield()`; irreversible
  -mortality + save-draws options), `make_temperature_scenarios()` (flat /
  single-spike / multi-spike / diurnal), `planted_dose_from_trace()` (analytical HI
  truth for validation).
- Plotting: `plot_survival_curves()`, `plot_tdt_curve()`, `plot_tdt_landscape()`,
  `plot_temperature_density()`, `plot_temperature_scenarios()`,
  `plot_heat_injury()`, `plot_repair_rate()`, shared `theme_tdt()`.
- Accessors / methods: `get_brmsfit()`, `has_fit()`, `get_z_draws()`,
  `get_ctmax_draws()`, `get_tcrit_draws()`, `get_surv_draws()`, `get_hi_draws()`,
  `print`/`summary`/`plot.bayes_tls`, `tdt_parameter_table()`,
  `summarise_observed_survival()`.

Files: `standardize_data.R`, `priors.R`, `fit_4pl.R`, `extract_tdt.R`,
`predict_survival_curves.R`, `predict_heat_injury.R`, `repair.R`,
`temperature_scenarios.R`, `tdt_landscape.R`, `plotting.R`, `accessors.R`,
`bayes_tls_methods.R`, `diagnostics.R`, `data.R`, `utils.R`,
`bayesTLS-package.R`. Tests mirror these in `tests/testthat/` (fast unit tests +
a `RUN_BRMS_TESTS=true` gated integration suite with recovery checks).

## The 4PL model

```text
p(t, T) = low + (up - low) / (1 + exp(k * (log10(t) - mid(T))))
```

with `low`/`up` the asymptotes (disjoint-bounds reparam above), `k = exp(logk) > 0`
the steepness on the log10-duration axis, and
`mid(T) = beta_0 + beta_1 * (T - T_bar) + …` the temperature-dependent midpoint.
Solving for duration at target survival `x`:

```text
log10(LTx(T)) = beta_0 + beta_1*(T - T_bar) + (1/k) * log((up - x)/(x - low))
```

For the fixed-shape, single-group case `z = -1/beta_1` exactly and the asymmetry
term shifts only the intercept. When shape parameters vary with `T`, the asymmetry
term bends the LT curve, so `z` becomes a local quantity (`extract_tdt()` handles
both via the relative/absolute threshold modes).

## Data

Five datasets ship as `.rda` in `data/`, built from raw CSVs in `inst/extdata/` by
`data-raw/make_datasets.R`, each documented in `man/`:

- `shrimp_lethal`, `shrimp_sublethal`, `zebrafish_lethal`, `snowgum_psii`,
  `acacia_seeds`.

`Rdata/` holds older shrimp intermediate objects from the legacy notebook; not
part of the current package-based workflow. Cached fits live in `output/models/`
(`fit_shrimp_lethal_4pl.rds`, `fit_zf_joint_4pl.rds`, `fit_leaf_function_4pl.rds`,
`fit_seed_lethal_4pl.rds`, the `sim_4pl_*` tutorial fits, etc.).

## Documents

- `ms/ms.qmd` — manuscript. Sections: introduction; the joint Bayesian 4PL model
  spec; deriving the TDT parameters; heat-injury/survival under dynamic
  environments; a simulation comparing the joint and two-stage approaches
  (`fig-sims`, `fig-sims2`; two-truth design); application to the case studies;
  discussion; conclusions. Conceptual figure is `fig-joint-4pl-conceptual`.
- `ms/supplement.qmd` — large (~5,500 lines): function tour, a simulated-data
  tutorial (two-stage pipeline vs joint 4PL, *z*/CTmax/*T*~crit~ derivation,
  T-varying-shape worked example, heat injury + planted-dose validation), the
  extended simulation results + sensitivity sweeps, and the four case-study
  sections. "Manuscript Figure 4" (cross-case-study summary) is generated here.

## Key decisions to remember (see CLAUDE.md memory + notes/)

- **Identity link** everywhere (counts and proportions) — validated equivalent to
  logit; safe because the disjoint-bounds reparam keeps the 4PL mean inside (0,1).
- **Beta family** is first-class for continuous proportions (leaf PSII). ~19% of
  leaf Fv/Fm are exact zeros, clamped to 0.001 — documented limitation.
- **`extract_tdt()` is not two-stage** — it is a per-draw posterior transform; the
  old per-draw OLS was removed (2026-05-24), `derive_tdt_parameters()` deleted,
  `derive_z()` added. Keep *z* and CTmax on the same `pars` subsample if you touch
  this path.
- **`temp_effects = "mid"`** for sparse designs — the all-four default can
  produce divergences and an unidentified *z* when cells, replication, or
  coverage of the survival range are limited.
- **Relative vs absolute survival threshold** — relative (mid-based) is the
  default; absolute is a sensitivity option (zebrafish §). They differ only when
  shape parameters vary with temperature.
- **Audience** is functional ecology / evolution / thermal physiology. Keep
  cross-disciplinary parallels (e.g. clinical CEM43) to a one-line citation; never
  reframe around them.

## Stale / superseded — ignore unless explicitly revived

- The old prototype notebooks are now in `notes/legacy/`: `bayesian_4pl.R`
  (zebrafish/life-stage prototype, expects a zebrafish xlsx not in the tree),
  `bayesian_TDT_shrimp.qmd` (broad shrimp notebook, several rough/session-
  dependent chunks), `test_functions.qmd` (the former reusable scaffold — its
  functions have been extracted, cleaned, and tested into `R/`). Use the package
  functions in `R/`, not these.
- Any reference to a Quarto book, `_quarto.yml`, `_book/`, `index.qmd`,
  `suppfig-*`/`supptbl-*` floats, or `R/bayesian_4pl.R`/`R/test_functions.qmd` as
  live code is from the retired architecture.

## Things not to forget

- `z = -1/b_mid_temp_c` is exact only for the fixed-shape single-group case; use
  the posterior-transform path for grouped/general models and non-50% thresholds.
- Repair parameters in `predict_heat_injury(repair = TRUE)` are scenario
  assumptions (Sharpe-Schoolfield defaults), not fitted from the TDT data.
- CTmax-like reference temperatures can be near-edge / extrapolated quantities;
  credible intervals expose the uncertainty but do not make extrapolation safe.
- Durations in the shrimp data are recorded in particular units; `standardize_data`
  records the unit and `extract_tdt`/`derive_tdt_curve` convert via
  `time_multiplier` — keep `t_ref` and the time axis consistent.
- Do not overwrite raw data; do not revert the user's modified/generated files
  without explicit permission.
</content>
