# Plan: layered direct-CTmax/z parameterisation in the package functions

**Status:** design / for discussion (precedes implementation).
**Date:** 2026-06-22.
**Builds on:** the verified maths in
`notes/2026-06-16-direct-ctmax-z-parameterisation.qmd` (relative + absolute
reparam, exactness vs the midpoint fit) and the prior-sensitivity finding in
`notes/2026-06-21-prior-sensitivity-direct-ctmax-z.qmd` (centred asymptote priors
are load-bearing).

## Goal (Daniel's framing)

Let users model **CTmax and z directly** (as the estimated parameters), while
**keeping** the ability to put fixed temperature effects and random effects on the
4PL shape parameters (`low`, `up`, `k`). A **layered** API: supplying CTmax/z
formulas switches on the direct parameterisation; a **clear `threshold` argument**
says whether CTmax/z are the *relative* or *absolute* quantities. The midpoint
parameterisation stays the default (backward compatible).

## Current architecture (what we're extending)

- `make_4pl_formula()` — 4PL is `low + (up-low)/(1+exp(exp(logk)*(logd - mid)))`.
  `low/up` are `inv_logit` reparams, `k = exp(logk)`, all inlined. The **estimated
  nlpars are `lowraw, upraw, logk, mid`**. `mid` is the linear-in-T backbone
  (`mid ~ temp_c + (1|re)`); `low/up/k` optionally get `~ temp_c` via
  `temp_effects`. Group moderators (life stage, sex) are *not* supported by
  `fit_4pl()` today — grouped fits are hand-built `bf()` (as in the case studies).
- `make_4pl_priors()` — centres `lowraw/upraw` so `low≈0.02`, `up≈0.98` (must keep
  — see prior-sensitivity note), `logk≈log2`, and a linear prior on `mid`.
- `fit_4pl()` — wraps formula+priors+`brm`, stores `meta` (temp_mean,
  duration_unit, bounds, temp_effects, family, ...).
- `tls()` / `extract_tdt()` — derive z, CTmax, T_crit from the **`mid` nlpar**:
  evaluate the four sub-parameters on a temperature grid via
  `posterior_linpred(nlpar=)`, take a per-draw LS slope of `log10 LT(T)`,
  `z = -1/slope`, `CTmax = temp_mean + (log10_tref - intercept)/slope`.

## The reparameterisation (from the verified note)

One linear-in-T backbone. `mid` is reconstructed from CTmax and z:

- **relative:** `mid(T) = log10_tref − (T_c − CTmaxdev)/exp(logz)`
  where `CTmaxdev = CTmax − temp_mean`, `z = exp(logz)`, `log10_tref` a fit-time
  constant (0 when tref = the model's time unit, e.g. 1 h in an hours model).
- **absolute (p-survival):** subtract the asymmetry correction
  `C(T) = (1/k)·log((up−p)/(p−low))` from the relative `mid`.

So the **estimated nlpars become `lowraw, upraw, logk, CTmaxdev, logz`** and `mid`
is a *derived* quantity. `z = exp(b_logz)`, `CTmax = temp_mean + b_CTmaxdev` are
read straight off the coefficients (no slope step) for the relative threshold.

## Proposed API (layered)

```r
fit_4pl(
  data,
  # --- existing shape-parameter layer (unchanged) ---
  temp_effects   = c("low","up","k","mid"),   # 'mid' meaning kept for midpoint param
  random_effects = NULL,                        # still attach to the backbone
  lower = 0, upper = 1, family = NULL, prior = NULL, ...,
  # --- NEW: direct CTmax/z layer ---
  parameterization = c("midpoint","direct"),   # default midpoint (back-compat)
  ctmax_formula    = NULL,                       # e.g. ~ life_stage + (1|Date)
  z_formula        = NULL,                       # e.g. ~ life_stage
  threshold        = c("relative","absolute"),   # meaning of CTmax/z when direct
  t_ref            = NULL                         # reference exposure (model units); default = 1 time-unit
)
```

**Trigger / layering.** Supplying `ctmax_formula` or `z_formula` ⇒
`parameterization = "direct"` (or set it explicitly). In direct mode:
- `CTmaxdev ~ <ctmax_formula>` and `logz ~ <z_formula>` carry the fixed +
  random structure for tolerance and sensitivity (this *replaces* the role
  `mid ~ ...` played; `random_effects`/group moderators now live here).
- `low/up/k` keep the **existing** `temp_effects` mechanism (and could later take
  their own formulas too) — that is the "keep effects on the 4PL parameters" half.
- `threshold` selects whether the linear backbone (and hence CTmax/z) is the
  relative midpoint or the absolute p-survival LT.

## Function-by-function change list

1. **`R/fit_4pl.R :: make_4pl_formula()`** — add `parameterization`, `threshold`,
   `ctmax_formula`, `z_formula`, `log10_tref`. When `direct`:
   - build `mid_expr` from CTmaxdev/logz (relative) or with `−C(T)` (absolute),
     inlined into the main 4PL (as the verified note does) **or** via
     `brms::nlf(mid ~ ...)` (cleaner; *verify `posterior_linpred(nlpar="mid")`
     still works on an `nlf`-derived param* — open question O1).
   - emit nlpar formulas `CTmaxdev ~ rhs(ctmax_formula)`,
     `logz ~ rhs(z_formula)`, keep `lowraw/upraw/logk ~ ...` as today.
   - `mid` is no longer an estimated nlpar.
2. **`R/fit_4pl.R :: fit_4pl()`** — new args; auto-detect direct from the
   formulas; compute `log10_tref` from `t_ref` + `meta$duration_unit`; record
   `parameterization`, `threshold`, `t_ref`, `ctmax_formula`, `z_formula` in
   `meta` (extractors branch on this). Validate: `threshold="absolute"` needs
   `low < p < up` to be attainable (guard/warn — O2).
3. **`R/priors.R :: make_4pl_priors()`** — add a `parameterization`/`threshold`
   branch. Direct priors: `CTmaxdev ~ normal(0, sd_ctmax)` (°C scale; default sd
   from the observed temp range, e.g. 10), `logz ~ normal(log(z0), 0.5–0.7)`
   (z0≈3). **Keep** the centred `lowraw/upraw/logk` priors unchanged (load-bearing
   per the prior-sensitivity note). Handle arbitrary `ctmax_formula`/`z_formula`
   coefficients with a single `class="b", nlpar=` prior (CTmax/z are
   well-identified, so generic-per-coef is safe here — unlike the asymptotes).
> **VERIFIED by the 2026-06-22 helper audit** (see §9 of
> `notes/2026-06-22-direct-ctmax-z-prototype-validation.qmd`): O1 is resolved —
> `posterior_linpred(nlpar="mid")` reaches the `nlf`-derived `mid`. So `tls()` +
> `tls_z/ctmax/tcrit`, `predict_survival_curves()`, `derive_tdt_landscape()`,
> `diagnose_tdt_fit()` and the accessors **work unchanged**. Only **three**
> functions break and need rewiring onto `posterior_linpred`: `extract_tdt()`
> (silently returns NA), `predict_heat_injury()` (hard error), and
> `tdt_parameter_table()` (hard error). Items 4–5 below are updated accordingly.

4. **`R/tls.R :: tls()`** — **no change required** (audit: correct values, the
   grid-slope on the `nlf`-derived `mid` equals the direct coef read). *Optional*
   ergonomic fast-path: when the fit is direct, read `z = exp(b_logz)` and
   `CTmax = temp_mean + b_CTmaxdev` directly instead of the slope — same answer,
   marginally cheaper.
5. **`R/extract_tdt.R :: extract_tdt()`** — **must change** (audit: silently
   returns all NA — it parses `b_mid_temp_c`/`b_mid_Intercept`). Rewire onto the
   `posterior_linpred` path, or have it delegate to `tls()`; for direct+relative,
   z/CTmax are coefficient reads; T_crit/absolute stay derived (the correction
   layer never disappears). Likely thin it to delegate to `tls()` where possible.
   Keep the existing `target_surv` arg as the **threshold converter**: it returns
   the fitted threshold (coef read) *or* the opposite threshold (via `C(T)`, i.e.
   `tls(mode = <other>)`) — validated both directions (prototype note §8/§9). So a
   relative fit can report absolute CTmax/z and vice versa; default to the fitted
   threshold (the linear backbone) and flag the converted one as derived.
8a. **`R/diagnostics.R :: tdt_parameter_table()`** — **must change** (audit: hard
   error, z column size 0; parses `b_mid_temp_c`). Same `posterior_linpred` fix.
6. **`R/bayes_tls_methods.R`** (`print`/`summary`) — report parameterisation,
   threshold, t_ref.
7. **`R/standardize_data.R`** — probably unchanged (t_ref is a fit arg); confirm
   `temp_mean`/`duration_unit` in meta cover what direct mode needs (they do).
8. **`R/diagnostics.R`** — `diagnose_tdt_fit()` lists nlpars; make it
   parameterisation-aware (CTmaxdev/logz vs mid).
9. **Tests (`tests/testthat/`)** — (a) direct≈midpoint equivalence on a shared
   dataset (z/CTmax medians + CIs match, as the note showed on zebrafish);
   (b) relative==absolute under constant shape (`temp_effects="mid"`-equivalent);
   (c) prior object well-formed for user formulas; (d) extractor agreement
   (coefficient read == slope derivation for relative).
10. **Docs / vignette** — roxygen for the new args + a worked example; fold the
    2026-06-16 tutorial note into the package tutorial when that lands.

## Key design decisions (need Daniel's call)

- **D1 — what `threshold` controls.**
  (A) it changes the **fitted backbone** (absolute bakes `−C(T)` into `mid`, so
  CTmax/z *are* the absolute quantities) — matches the note's `mid_absolute`;
  vs (B) always fit the **relative** backbone and only *derive* absolute at
  extraction (current philosophy; simpler, better-conditioned). Recommendation:
  default to modelling **relative** directly; offer `threshold="absolute"` as
  (A) for users who want absolute CTmax/z as parameters, with the well-posedness
  guard. Absolute LTx/T_crit remain derivable either way.
- **D2 — trigger.** Auto-switch to direct when `ctmax_formula`/`z_formula` are
  supplied (ergonomic) vs require explicit `parameterization="direct"` (explicit).
  Recommendation: accept the formulas as the trigger **and** allow the explicit
  arg; error if they conflict.
- **D3 — `mid` via `nlf()` vs inlined.** `nlf()` is cleaner and may keep
  `posterior_linpred(nlpar="mid")` working (so the derivation-based extractors and
  absolute/T_crit need little change); inlining is what the note verified. Decide
  after testing O1.
- **D4 — effects on low/up/k.** Keep the simple `temp_effects` toggle now, or also
  let users pass `low_formula`/`up_formula`/`k_formula` for full symmetry?
  Recommendation: ship `ctmax_formula`/`z_formula` first; add shape formulas later
  if needed (keeps scope contained).
- **D5 — default `t_ref`.** 1 hour (matches the manuscript's CTmax₁ₕ) expressed in
  the model's time unit, i.e. `log10_tref = log10(1h / duration_unit)`.

## Open implementation questions

- **O1** Does `posterior_linpred(fit, nlpar="mid")` work when `mid` is an
  `nlf()`-derived (not estimated) parameter? If yes, extractors stay nearly
  unchanged; if no, derive the backbone from `CTmaxdev/logz/low/up/k` draws.
- **O2** Absolute threshold well-posedness when `up < p` (sublethal) — warn and
  fall back to relative, or error.
- **O3** Sampling geometry of `/exp(logz)` (funnel when z weakly identified) on
  sparse designs — carry the note's caveat; test ESS/divergences.

## Backward compatibility

`parameterization="midpoint"` (default) leaves `make_4pl_formula`,
`make_4pl_priors`, `fit_4pl`, `tls`, `extract_tdt` behaving exactly as now. Direct
mode is purely additive. Ship as a minor/major version bump; keep both
parameterisations indefinitely (the correction layer means the midpoint form is
still the natural home for absolute/LTx/T_crit).

## Validation gate before merge

Reproduce the note's zebrafish equivalence through the *package functions* (not a
hand-built `bf`): `fit_4pl(parameterization="direct", ctmax_formula = ~ life_stage
+ (1|Date_experiment), z_formula = ~ life_stage + (1|Date_experiment))` must match
the midpoint joint fit's z/CTmax (medians ~0.01°C; CIs match when RE structures
match), with 0 divergences under the centred asymptote priors.
