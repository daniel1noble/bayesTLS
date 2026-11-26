# bayesTLS (development version)

## New features

* Censored continuous-proportion readings. `standardize_data()` gains
  `detection_limit` for a `proportion` response whose instrument cannot
  resolve values beyond a limit (e.g. a fluorometer that records 0 below
  what it can read, or values capped at a ceiling). Give it as one number
  (a lower limit), `c(lower, upper)`, or per-row limit columns. Readings
  below the lower limit are flagged left-censored and readings at or above
  the upper limit right-censored (a new `censored` column, with the limit as
  the response). The fitters (`fit_4pl()`, `fit_lognormal()`,
  `fit_weibull()`) then add brms's `cens()` term to the Beta response, so
  those readings contribute the probability of lying beyond the limit
  instead of a density at a clamped value; `make_4pl_formula()` and
  `make_tdt_formula()` gain `censored`. The fitted curve keeps its meaning,
  so `tls()` and every prediction helper are unchanged.
* Plots (`plot_survival_curves()`, `plot_tdt_landscape()`) and
  `summarise_observed_survival()` show censored readings as recorded
  (`survival_recorded`), `print()` reports the censored counts, and
  `bayes_R2_tls()` warns that R^2 is only a rough guide when readings are
  censored. A censored fit scores censored readings as probabilities in
  `log_lik()`, `loo()` and `waic()`, so compare it only with fits to the same
  data and the same `detection_limit`, never with a fit whose readings were
  clamped. The new arguments come last in each signature, so existing
  positional calls are unchanged. Without `detection_limit` nothing changes,
  except that `proportion_eps` must now be one number in (0, 0.5) (values
  outside that range gave degenerate or invalid responses before).

## Bug fixes

* `derive_tdt_curve()` with an absolute threshold now reads each draw's LT
  exactly from its fitted curve: the midpoint plus the model's closed-form
  correction, the same values `tls()` and `derive_z()` use. It used to
  search a duration grid and silently drop every draw whose crossing lay
  beyond it, so cool temperatures were summarised from a biased few draws
  (on the snow gum fit, dark recovery at 40 degrees C came from 1 of 400
  draws). `duration_grid` is no longer used, and a message says so if it is
  given. Draws whose curve never reaches the target (outside its
  asymptotes) are dropped with a warning, as elsewhere.
* `tls()` and `derive_temperature_for_duration()` no longer silently drop
  draws whose CTmax lies outside the temperature grid, which made intervals
  too narrow. When the LT curve is a straight line in temperature, CTmax
  (and the 1-hour anchor of T_crit) now come from the exact closed form at
  any temperature. That covers every relative threshold, for midpoint and
  direct fits alike, and absolute thresholds whose asymptotes and steepness
  do not change with temperature. `tls()` used to invert direct fits on the
  grid, dropping draws even at the default reference exposure. A bent curve
  (an absolute threshold with temperature effects on the asymptotes or
  steepness) is still inverted on the grid, and draws that do not cross
  there are now searched outward, up to the data's temperature range (or
  the grid's own width, if larger) beyond each end, taking the crossing
  nearest the data (an extrapolated bent curve can turn back and cross
  again further out). Draws still without a crossing are counted, per
  quantity, in a warning. The grid is first refined so that neighbouring
  points are at most 0.05 degrees C apart, or a twentieth of the original
  spacing where that is wider than 1 degree C (the snow gum default's
  2.6-degree steps become 0.13): interpolating linearly between 11 coarse
  points could put a crossing up to about 2 degrees C from the curve's own
  crossing next to an asymptote. Against a brute-force search on a fine
  grid, the largest per-draw error fell from 0.26 to 0.03 degrees C on the
  test fits and from 2.1 to 0.06 degrees C on a censored snow gum fit, with
  medians within 0.002 degrees C. In the supplement, the zebrafish normoxia
  CTmax moves by 0.01 degrees C (10 of its 8000 draws had been dropped) and
  every other number by less than 0.003 degrees C. The search evaluates a
  group's draws in bounded chunks and keeps, between its two passes, only
  the draws still without a crossing, so its memory no longer grows with
  the number of groups.
  `derive_temperature_for_duration()` no longer assumes
  a straight relative midpoint (one that bends is inverted instead), and
  its `temp_grid` now defaults to 11 points over the observed temperatures.
  Its check for a straight curve no longer looks only at the requested
  temperatures, so a `temp_grid` of the centring temperature and/or one
  degree above it cannot make a bent curve look straight (it used to be
  solved as the line through those two temperatures, off by up to 2 degrees
  C per draw on a test fit). A one-point `temp_grid` is searched outward at
  a fine step instead of in quarters of the data's temperature range, and a
  grouped call in which one group has no crossing gives that group an `NA`
  summary row, with the warning, instead of an error.
* Absolute thresholds near an asymptote (for example `target_surv = 0.9`
  when the upper asymptote falls with temperature): where the target lies
  outside a draw's fitted asymptotes at some temperatures, the LT there is
  undefined, and the CTmax inversion used to drop every draw with such a
  temperature on its grid. It now treats that LT as zero (upper asymptote
  at or below the target) or infinite (lower asymptote at or above it) and
  locates a crossing next to it from the survival curve, which is defined
  everywhere. On a test fit at `target_surv = 0.8`, 1000 of 1000 draws now
  have a CTmax, against 2 before.
* The grid inversion of a bent LT curve that is not monotone over the grid
  (U- or cap-shaped near an asymptote) used to interpolate temperature
  against the sorted LT values, which mixes the two branches and can return
  a temperature that is not a crossing at all (by up to about 10 degrees C
  for some draws of a real fit at `target_surv = 0.9`). It now takes the
  single crossing where the LT shortens through the target as temperature
  rises; a draw with none, or more than one, gets no CTmax and is counted
  in the warning. Relative thresholds and `target_surv = 0.5` were not
  affected on any fit checked.
* `tls()` with `ndraws` now uses one posterior subsample for all four 4PL
  sub-parameters. Each of the four (low, up, k, mid) used to come from its
  own random subsample, so they were mismatched draws, which distorted
  absolute-threshold results (relative thresholds use the midpoint only and
  were unaffected). Seeded results with `ndraws` change, since the
  subsample is drawn differently. An `ndraws` at least the posterior size
  now uses all draws instead of erroring, and every reader with an
  `ndraws` argument (including `predict_survival_curves()`,
  `derive_tdt_landscape()` and `predict_heat_injury()`) now requires it to
  be `NULL` or a whole number of at least 1.
* `tls()` warns when `temp_grid` looks like degrees C rather than centred
  temperatures: when it lies entirely outside the observed (centred)
  temperatures, or reaches more than the data's own temperature range
  beyond them. A degrees-C grid used to give nonsense silently. The check
  cannot tell the two scales apart for data near 0 degrees C. Non-finite
  `temp_grid` values are now dropped with a warning in `tls()`,
  `derive_tdt_curve()` and `derive_temperature_for_duration()`; an `Inf`
  used to make the crossing search fail silently.
* The warning for an absolute `target_surv` outside the fitted asymptotes
  now gives the number of draws affected and on which side of the curve
  the target falls. It used to round the share of evaluations (it could
  say "0%" while dropping some) and print asymptote medians over all of
  them.
* `derive_tdt_curve()`, `derive_z()`, `derive_temperature_for_duration()`,
  `tdt_parameter_table()`, `extract_4pl_pars()` and the heat-injury
  functions now work on a hand-built `bayes_tls` that lacks `meta$bounds`,
  rebuilding the bounds from `meta$lower` and `meta$upper` as `tls()` does.
  Absolute thresholds used to return nothing, with a misleading warning.
* `tls()` with `by = NULL` on a grouped fit (e.g. `ctmax = ~ 0 + recovery`)
  now reports every group, as `derive_z()` and the other readers already
  did. Previously it filled the moderator from the first data row and
  silently returned that one group without a group label. A bare `brmsfit`
  whose CTmax/z (or midpoint) vary by a moderator now asks for `by`; with a
  user-supplied `newdata` grid, `by = NULL` still treats the grid as one
  group. Single-condition fits are unaffected.

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
  fixed survival target), `predict_heat_injury()` (dose and survival from
  the fitted model's curve: read off it for the constant shape, followed per
  individual or by cumulative exposure for the varying shape), survival curves,
  landscapes and `tdt_parameter_table()` (which records the model in its
  `"model"` attribute). A bare `brmsfit` passed to `tls()` is recognised from
  its formula.
* Log-logistic results are unchanged. On the paper's 12 cached fits, `tls()`
  z and CTmax (relative and absolute thresholds) and T_crit, as summaries
  and as every posterior draw, the curve parameters, `derive_z()`,
  `derive_temperature_for_duration()` and `derive_tdt_curve()` at the
  absolute threshold, and `predict_heat_injury()` with the default constant
  shape (relative and absolute dose) are bit-identical to 1.1.0, compared
  on the centred temperature scale and again with a shifted centring
  constant; `tdt_parameter_table()` gains only the `"model"` attribute
  (`scripts/package/check_case_study_outputs.R`). The varying shape changes:
  see Bug fixes.
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

* `predict_heat_injury(shape = "varying")` used an ad hoc survival update,
  present since the shape was added in 1.0.0. The default
  `shape = "constant"`, which the paper uses, gives the same numbers for
  every draw that produces a prediction (draws that do not are now blanked;
  see the next entry). Under the old update:
  - each step subtracted the fixed number of points the current curve drops,
    whatever share of the animals was still alive, so survival could fall
    below every curve's lower asymptote (and was clipped at 0), could rise
    when time at a hot temperature was added, and depended on the time of day
    a trace started;
  - with an absolute `target_surv`, steps where the target was at or above
    the curve's upper asymptote added no heat injury and no mortality, and a
    draw whose first step was such a step returned `NA` and was dropped from
    the summary without a warning;
  - with repair and irreversible mortality, injury that was repaired and then
    re-accrued was charged again (identical repaired days kept killing). That
    is what stochastic death predicts, and it was inconsistent with the
    individual-tolerance rule of the constant shape, under which injury kills
    no one until it passes its previous peak.

  The varying shape now offers two coherent methods, chosen with the new
  argument `varying_method`: `"individual"` (default), where each individual
  accumulates heat injury additively against its own threshold time, and
  `"cumulative"`, cumulative exposure with survival as the state (the rule
  used by Rezende et al. 2020 and Bullard et al. 2026; given the same curves
  it reproduces Bullard et al.'s code up to the time step). Both reproduce the
  fitted curve at any constant temperature and both equal the constant shape
  exactly when `low`, `up` and `k` do not change with temperature and `mid`
  is linear in it, with or without repair; they differ otherwise, because
  constant-temperature data cannot say why the shape changes. Heat injury is
  the same under both unless repair scales with survival (for an absolute
  target, the injury of the individual of that rank; `Inf`, with a warning,
  while injury accrued where the target is already passed is outstanding).
  Both read a lower `up` at a temperature as individuals that die on exposure
  there, permanently even with reversible mortality; under cumulative
  exposure a temperature whose floor is at or above survival adds no injury. Survival starts at `up` of the first temperature
  that accrues injury (above `T_c`; a trace that accrues none keeps the
  constant shape's value), and `repair_scales_with_survival`
  scales repair by survival over that `up`, a factor never above 1. The help
  page sets
  out both methods and their assumptions, and the shape is marked
  experimental. Varying-shape results change for fits whose `low`, `up` or
  `k` change with temperature, for any fit with repair and irreversible
  mortality, and for absolute targets outside a draw's curve at some trace
  temperature.
* `predict_heat_injury()` now reports in `meta$ndraws` the posterior draws
  that produced a prediction, and warns when any did not (for example, an
  absolute target outside a draw's curve at the centring temperature). Such
  draws are now left out of every summary column, and `draws` keeps them as
  rows of `NA`: before, the count reported was the number requested, the
  draws were dropped silently, and with `repair_scales_with_survival = FALSE`
  their heat injury still entered the summary as 0. For a grouped fit,
  `meta$ndraws` is the smallest count in any group and the new
  `meta$ndraws_by_group` gives each group's count.
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
