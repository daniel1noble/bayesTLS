# bayesTLS 1.0.0

First release.

`bayesTLS` fits joint Bayesian four-parameter logistic (4PL) models to
thermal-tolerance data and derives the classical thermal death time / thermal
load sensitivity quantities from the posterior, so that every downstream
quantity carries full uncertainty and is mutually consistent within a draw.

## Model fitting

* `standardize_data()` maps a raw thermal-tolerance dataset (binomial counts or
  continuous proportions) onto the columns the model expects, and records the
  duration unit and centring used.
* `fit_4pl()` fits the joint 4PL with `brms`, in either the *midpoint*
  parameterisation or the *direct* CTmax/z parameterisation
  (`ctmax = ~ ..., z = ~ ...`), with moderators allowed on any sub-parameter.
* `make_4pl_formula()` and `make_4pl_priors()` expose the underlying `brms`
  formula and default priors for inspection or customisation.

## Deriving thermal load sensitivity

* `tls()` derives `z`, `CTmax` and `T_crit` per moderator group from any fitted
  4PL — including hand-written `brms` models — by evaluating the sub-parameters
  on a moderator x temperature grid with `brms::posterior_linpred()`.
* `tls_z()`, `tls_ctmax()` and `tls_tcrit()` return the individual quantities;
  `derive_tdt_curve()` and `derive_tdt_landscape()` give the TDT curve and
  landscape.

## Prediction under fluctuating temperatures

* `predict_heat_injury()` accumulates heat injury over an arbitrary temperature
  trace, with optional Sharpe-Schoolfield repair (`repair_rate_schoolfield()`).
* `predict_survival_curves()` propagates that injury to survival.
* `make_temperature_scenarios()` builds fluctuating-temperature scenarios.

## Two-stage comparison

* `ts_stage1()`, `ts_stage2()`, `ts_curve()` and `ts_ci()` implement the
  conventional two-stage TDT workflow, for direct comparison against the joint
  model.

## Plotting and diagnostics

* `plot_tdt_curve()`, `plot_tdt_landscape()`, `plot_heat_injury()`,
  `plot_survival_curves()`, `plot_temperature_scenarios()` and friends, all on a
  shared `theme_tdt()`.
* `diagnose_tdt_fit()` and `bayes_R2_tls()` for fit checking.

## Data

Seven example datasets spanning lethal and sub-lethal endpoints: `aphid_tdt`,
`dsuzukii`, `shrimp_lethal`, `shrimp_sublethal`, `snowgum_psii`,
`zebrafish_lethal` and `zebrafish_o2`.
