# bayesTLS function map audit

Date: 2026-06-25

Artifacts:

- `scripts/make_function_map.R`
- `output/figs/bayesTLS_function_map.svg`
- `output/figs/bayesTLS_function_map.png`

## Scope

The map is deliberately a core-use map, not a full call graph. I inventoried the
exported surface from `NAMESPACE`, top-level functions in `R/`, tests, and the
function tour in `ms/supplement.qmd`. The full package has many exported helpers,
but the user-facing workflow can be explained with a much smaller spine:

1. `standardize_data()`
2. `fit_4pl()`
3. `tls()` and/or `extract_tdt()`
4. `predict_survival_curves()` and `predict_heat_injury()`
5. `get_tls_est()`, `get_4pl_est()`, and the `plot_*()` family

Everything else is either advanced model specification, a primitive used by a
wrapper, diagnostics, plotting support, a classical two-stage comparison path, or
an accessor/convenience wrapper.

## Recommended documentation tiers

**Core route.** Use these in the first-page map and the main tutorial:
`standardize_data()`, `fit_4pl()`, `tls()`, `predict_survival_curves()`,
`predict_heat_injury()`, `get_tls_est()`, `get_4pl_est()`, and grouped `plot_*()`
families.

**Workflow-specific route.** Keep `extract_tdt()` visible, but explain it as the
`fit_4pl()` workflow bundle that returns the historical nested object plus an
LT curve. It is useful, but it should not compete with `tls()` as the general
mental model.

**Advanced route.** Put `derive_z()`, `derive_tdt_curve()`,
`derive_temperature_for_duration()`, `make_4pl_formula()`, `make_4pl_priors()`,
`tdt_parameter_table()`, and `extract_4pl_pars()` in an advanced/API reference
section. They are real functions, but they are not first-pass workflow verbs.

**Comparison route.** Keep `ts_stage1()`, `ts_stage2()`, `ts_ci()`, and
`ts_curve()` as a separate two-stage comparator module. They should not sit in
the same visual lane as the Bayesian workflow.

## Consolidation candidates

### 1. Accessor proliferation

Current overlap:

- Individual `extract_tdt()` accessors: `get_z_draws()`, `get_z_summary()`,
  `get_ctmax_draws()`, `get_ctmax_summary()`, `get_tcrit_draws()`,
  `get_tcrit_summary()`
- Combined `extract_tdt()` accessors: `get_tls_draws()`, `get_tls_summary()`
- `tls()` accessor: `get_tls_est()`

Recommendation: make `get_tls_est()` the main public accessor for both `tls()`
and `extract_tdt()` objects. Keep the individual functions as back-compatible
wrappers, but move them out of the core docs and consider soft-deprecation once
the manuscript/tutorial code stops using them directly.

Reason: the current interface teaches users two output ecosystems for the same
scientific quantities. The combined accessor already expresses the better mental
model: one tidy object of quantities, summaries, or draws.

### 2. `tls()` versus `extract_tdt()`

`tls()` is the cleaner general verb: it works on a `bayes_tls` workflow or a
hand-coded `brmsfit`, returns tidy-long summaries/draws, and supports
moderator-group extraction. `extract_tdt()` remains useful because it preserves
the established nested contract, returns the descriptive LT curve, and performs
the workflow-specific true inversion/local-z path.

Recommendation: do not merge them immediately. Promote `tls()` as the primary
extraction verb in docs, and describe `extract_tdt()` as a workflow-specific
bundle/compatibility function. A later merge is only safe after the absolute
threshold semantics are reconciled: `tls()` currently reports a least-squares
slope summary over the LT curve, while `extract_tdt()` uses local z and true
inversion.

### 3. Derivation primitives

`derive_z()`, `derive_tdt_curve()`, and `derive_temperature_for_duration()` are
useful primitives, but they duplicate conceptual territory with `extract_tdt()`
and `tls()`.

Recommendation: keep them exported for advanced users, but keep them out of the
core workflow figure. If the API is simplified later, consider a single
`derive_tdt()`/`derive_threshold()` style function with `quantity = c("z",
"curve", "temperature")`, while retaining wrappers for compatibility.

### 4. 4PL parameter accessors

`get_4pl_est()` is already the general public wrapper over
`tdt_parameter_table()` and `extract_4pl_pars()`.

Recommendation: make `get_4pl_est()` the documented public route. Treat
`tdt_parameter_table()` and `extract_4pl_pars()` as advanced/internal-ish
building blocks unless there is a strong user-facing need for both names.

### 5. Convenience wrappers

`tls_z()`, `tls_ctmax()`, and `tls_tcrit()` are one-line wrappers around
`tls(params = ...)`.

Recommendation: omit them from the core map and the main tutorial. They can stay
for convenience/back-compat, but they add more names than concepts.

### 6. Prediction wrappers and validation helpers

`derive_tdt_landscape()` is a useful grid preset around
`predict_survival_curves()`. `summarise_observed_survival()` is mostly a plotting
overlay helper. `planted_dose_from_trace()` is a validation/tutorial truth
calculator.

Recommendation: keep `derive_tdt_landscape()` in plotting/prediction examples
where a heatmap is needed, but not in the first-route map. Consider whether
`summarise_observed_survival()` and `planted_dose_from_trace()` need to be
exported at all, or whether they belong in tests/tutorial support.

## Low-risk next cleanup

The least disruptive simplification would be documentation-first:

1. Keep all current exports for now.
2. Replace the full function table's prominence with the compact function map.
3. Teach `get_tls_est()`/`get_4pl_est()` as the preferred accessor layer.
4. Move individual `get_z_*`, `get_ctmax_*`, `get_tcrit_*`, `tls_*` wrappers,
   and validation helpers into an "advanced and compatibility helpers" section.

That reduces cognitive load without breaking existing manuscript code or tests.
