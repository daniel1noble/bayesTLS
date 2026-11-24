# bayesTLS 1.2.0

## New features

* Alternative dose-response models with the same four parameters and the
  same TDT machinery: `fit_lognormal()` (log-normal / probit curve) and
  `fit_weibull(type = 1 | 2)` (Weibull type 1, the Gompertz-shaped curve on
  log time, and Weibull type 2). Every model is written in its midpoint
  form, so `mid`, the TDT line, the relative threshold, `z` and `CTmax` mean
  the same thing whichever curve is fitted; the prior mean of `logk` is
  shifted per model so the prior on the slope at the midpoint matches.
  `fit_4pl()` is unchanged (it keeps meaning the log-logistic and gains no
  model argument); `make_tdt_formula()` / `make_tdt_priors()` are the
  model-generic builders behind `make_4pl_formula()` / `make_4pl_priors()`.
  Everything downstream works for every model: `tls()`, `derive_z()`,
  `derive_temperature_for_duration()` and `derive_tdt_curve()` at relative
  or absolute thresholds (each model's own correction from the midpoint to a
  fixed survival target), `predict_heat_injury()` (dose and survival read
  off the fitted model's curve, constant or varying shape), survival curves,
  landscapes and `tdt_parameter_table()` (which records the model in its
  `"model"` attribute). A bare `brmsfit` passed to `tls()` is recognised from
  its formula.
* Log-logistic results are unchanged. On the paper's 12 cached fits, `tls()`
  z and CTmax (relative and absolute thresholds) and T_crit, as summaries
  and as every posterior draw, the curve parameters, `derive_z()`,
  `derive_temperature_for_duration()` and `derive_tdt_curve()` at the
  absolute threshold, and `predict_heat_injury()` (constant and varying
  shape, relative and absolute dose) are bit-identical to 1.1.0, compared
  on the centred temperature scale and again with a shifted centring
  constant; `tdt_parameter_table()` gains only the `"model"` attribute
  (`scripts/package/check_case_study_outputs.R`).
* `fit_4pl(fix_asymptotes = )` (and `make_4pl_priors()`) can pin one or both
  asymptotes instead of estimating them: `TRUE` fixes `low` and `up` at the
  response bounds (0 and 1 for proportions, a two-parameter log-logistic),
  `c("up")` fixes only the upper asymptote while `low` stays estimated, and
  `c(low = 0.02)` pins a named asymptote at a chosen value. Pinning is done
  through a brms `constant()` prior on the asymptote's Intercept with its
  formula forced to `~ 1`, so `tls()`, `predict_heat_injury()` and the
  accessors read a fixed asymptote exactly as an estimated one. Values at the
  response bounds sit at the edge of the padded window (0.001 / 0.999 by
  default); the effective values are recorded in `meta$fix_asymptotes` and
  shown by `print()`. The default (`NULL`) is unchanged.

## Bug fixes

* `tls()` now reads a workflow's own asymptote bounds (the ones it was fitted
  with) when it computes absolute thresholds, as `derive_z()` and the other
  readers already did. Before, it always used the default `[0, 1]` bounds, so
  absolute-threshold summaries of a workflow fitted with other bounds (e.g.
  `lower = 0.85` for photosystem-II data) were wrong without warning;
  relative thresholds were never affected. Fits with the default bounds,
  including all of the paper's, give identical results. `lower` and `upper`
  still override, and are still used for a bare `brmsfit`.

# bayesTLS 1.1.0

## New features

* Covariance-matrix random effects (issue #12): `fit_4pl()` now accepts
  `brms`'s `gr()` grouping terms, e.g.
  `ctmax = ~ 1 + (1 | gr(species, cov = A))` with a phylogenetic correlation
  matrix supplied through the new `data2` argument (`ape::vcv(phy, corr =
  TRUE)` convention; validated for symmetry, level-covering dimnames and unit
  diagonal). Internally, random-effect terms are now parsed from the formula
  syntax tree instead of with text patterns, which previously broke on the
  nested parentheses of `gr()` calls.
* `tls()` auto-includes matching random-effect terms in
  `posterior_linpred()` when a `by` variable enters the model only as a
  random-effect grouping factor (with a message; explicit `re_formula`
  always wins). Previously such calls silently returned the population mean
  for every group. All existing fits -- where `by` is a fixed-effect
  moderator -- are unaffected.

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

Four publicly available example datasets spanning lethal and sub-lethal
endpoints: `aphid_tdt`, `dsuzukii`, `snowgum_psii` and `zebrafish_o2`.
