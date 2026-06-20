# Audit — `profileTLS` (itchyshin/profileTLS) vs `bayesTLS`

**Date:** 2026-06-16
**Auditor:** Claude (package-maintenance/statistician role)
**Target:** https://github.com/itchyshin/profileTLS @ `6f963a9` (v0.3.3), cloned read-only to `/tmp/profileTLS`.
**Brief:** full audit of code, calculations, statistics; errors and uncertainty; alignment with the bayesTLS models; are results the same; how much time saved.
**Scope of this audit (honest):** I read and verified the *statistically load-bearing* code myself — the TMB C++ likelihood (`src/profile_tls.cpp`), `heat_injury.R`, `fit_engine.R`, and the equivalence algebra. The inference machinery (`profile.R`, `confint.R`, `bootstrap.R`, `diagnostics.R`), the data-build/benchmark/timing/coverage artefacts, and the bayesTLS-manuscript data path were verified via delegated reviewers who cited file:line and pasted real R output (live re-fits + cache loads); I cross-checked their headline claims. I did **not** deeply review `formula.R`, `plotting.R`, `simulate.R`, `methods.R`, `model_matrix.R` (lower statistical risk). No files in either repo were edited.

---

## Headline verdict

`profileTLS` is a **careful, competent maximum-likelihood / profile-likelihood reimplementation** of the bayesTLS 4PL TLS model. The model and likelihood are **mathematically equivalent to bayesTLS** and the calculations check out. For well-identified count data it reproduces bayesTLS point estimates to **~0.06 °C (CTmax) and <1 % (z)**, roughly **5× faster** (more for the Wald path).

Three things matter most:

1. **It surfaced a genuine data bug in *bayesTLS itself* (the shrimp dataset) — and that bug is in our shipped data AND our manuscript.** This is the single most important finding and is independent of profileTLS's own quality.
2. **profileTLS has several real inference rough edges** (silent non-convergence; an `up` confidence interval that can fall outside [0,1]; beta-binomial profile CIs that under-cover).
3. **The actual code has grown well beyond its own SPEC's v0.1 non-goals** (Beta/continuous, random effects, heat injury are all implemented now), and the newer additions carry the weak spots.

---

## 1. CRITICAL — a real bug in **bayesTLS** (our shrimp data), surfaced by profileTLS

This is about our package, not profileTLS's. profileTLS found it, documented it ("R-SHRIMP"), and fixed it on their side.

- **The bug:** `data-raw/make_datasets.R:34` runs `Mortality_after_trial = as.integer(Mortality_after_trial)`. But in the source CSV (`inst/extdata/data_lethal_TDT_brown_shrimp.csv`) `Mortality_after_trial` is a **proportion** (e.g. `0.0909 = 1/11`, `0.5 = 5/10`), and the comment at `:25` mislabels it a "death count". `as.integer()` truncates every proportion `< 1` to `0`.
- **Shipped object:** `data/shrimp_lethal.rda` has `Mortality_after_trial` collapsed to `{0: 113 rows, 1: 35 rows}` → **35 deaths / 1499 individuals = 2.3 %**.
- **What the manuscript actually fits:** `ms/supplement.qmd:2152-2162` feeds the shipped (corrupt) column to `standardize_data()`, which does `n_surv = round((1-mortality)*n_total)`. Because the surviving 35 rows are `mortality==1` (whole-tank kills) and 113 are `0`, the response entering the model is a clustered all-or-nothing pseudo-dataset at **23.7 % mortality** (`n_dead ∈ {0,10,11}`). The cached fit `output/models/fit_shrimp_lethal_4pl.rds` was trained on exactly this (`n_dead` frac 0.2368). The rendered supplement reports shrimp `z = 2.35 [1.76, 3.26]`, `CTmax_1hr = 32.66 [32.22, 33.40]`, `T_crit = 26.85 [24.03, 28.69]` — **all on corrupted counts.**
- **Correct data:** rebuilding from the CSV as `deaths = round(prop * N)` gives **49.2 % mortality, deaths spanning 0–11** (738 deaths). profileTLS ships this corrected version and guards against re-collapse (`make_benchmark_data.R:32-47, 157-161`).

**Why it produced "plausible" numbers anyway:** 23.7 % all-or-nothing mortality still spans a usable dose-response, so the fit converged to sensible-looking `z`/`CTmax` — the two-stage and joint methods "agree" only because they agree *on the same corrupt input*. This is exactly the kind of bug that does not announce itself.

**Action (high priority, our side):** fix `make_datasets.R` (`deaths = round(prop * N)`), rebuild `shrimp_lethal.rda`, refit `fit_shrimp_lethal_4pl.rds`, re-render the shrimp case study, and check whether the headline shrimp `z`/`CTmax`/`T_crit` move. Credit Shinichi's team for catching it. (Note: this is also the kind of `data/` provenance check CLAUDE.md §9 calls for.)

**Consequence for the cross-package benchmark:** profileTLS's benchmark feeds the **corrected** shrimp data to `bayesTLS::fit_4pl()` (`build_benchmark_cache.R:80-96`). So the benchmark's bayesTLS column and profileTLS agree because *both used corrected data* — but that is **not the data our manuscript uses**. Any "profileTLS == bayesTLS on shrimp" statement is on corrected data; our published shrimp fit is on the corrupt 23.7 % data. They are not the same dataset.

---

## 2. Calculations & alignment with the Bayesian model — **verified correct**

- **Reparameterisation equivalence (algebra checked).** profileTLS parameterises the midpoint directly: `mid_i = log10(tref) − (temp_i − CT_i)/z_i` (`src/profile_tls.cpp:167`). Expanding in centred form gives slope `β1 = −1/z` and `mid = log10(tref)` exactly at `temp = CT`. This is **identical** to bayesTLS's `mid(T)=β0+β1(T−Tbar)`, `z = −1/β1`, `CTmax(tref)=Tbar+(log10(tref)−β0)/β1`. Same `(low,up,k)`; smooth invertible reparam ⇒ same likelihood, curve, and MLE — but CTmax and z are now coordinates, hence directly profile-able. This is the legitimate, elegant core idea.
- **4PL and likelihoods (C++ read line-by-line).** `p = low + (up−low)·invlogit(−k(logd−mid))` = bayesTLS's `low + (up−low)/(1+exp(k(logd−mid)))` ✓ (descending in duration). Binomial `dbinom` on survivors ✓. Beta-binomial via the lgamma form with `a=p·phi, b=(1−p)·phi` is the standard PMF with `phi = a+b` precision — **matches brms's `beta_binomial` convention** (R-PHI is genuinely controlled) ✓. Beta family uses the same precision convention as brms `Beta(identity)` ✓. Probability clamp and shape floors are sound.
- **Heat-injury integral aligns with bayesTLS.** `LT(T) = tref·10^((CTmax−T)/z − q/k)`, `dmg = 1/LT`; at the relative midpoint (`q=0`, `tref=1`) `dose = ∫10^((T−CTmax)/z)dt`, `injury = 100·dose` — exactly our `HI = 100∫10^((T−CTmax)/z)di`. Forward-Euler with real per-step `Δt`, `T_c` cutoff (zero damage below threshold), monotone survival, and an optional Sharpe–Schoolfield repair kernel in Kelvin — all match `bayesTLS::predict_heat_injury()`. It's a faithful deterministic/ML analogue (point estimate + parametric-bootstrap band vs our posterior band).
- **One genuine modelling difference to note:** profileTLS uses a **nested-gap** asymptote reparam (`low=plogis(β_low)`, `up=low+(1−low)·plogis(β_gap)`), whereas bayesTLS uses **disjoint bounds** forcing `low<0.5<up`. Both guarantee `up>low`; profileTLS's is more flexible. On well-behaved data (low≈0, up≈1) they coincide, but on data where bayesTLS's bounds bind (asymptotes near 0.5) the two can give different `low`/`up` and therefore different absolute thresholds. Not a bug — a deliberate, defensible divergence — but it means "same model" is true only in the interior.

---

## 3. Are the results the same?

| Case | CTmax agreement | z agreement | Verdict |
|---|---|---|---|
| shrimp (lethal, binomial-ish) | 0.057 °C | 0.6 % | ✓ tight — **but on corrected data, not the manuscript's** |
| *D. suzukii* F/M (lethal) | 0.028 / 0.029 °C | 0.3 / 0.1 % | ✓ tight; reproduces Ørsted CTmax≈35.2, z≈3.0/3.2 |
| snowgum PSII (**Beta**) | **0.567 °C** | **~12 %** | ✗ outlier — exceeds any "tight agreement" bar |

- **Point estimates** match bayesTLS to high precision for well-identified **count** data (verified by live re-fit, not just cache). This is the expected consequence of the exact likelihood equivalence.
- **The snowgum PSII (continuous Beta) case is a real disagreement** (0.57 °C CTmax, ~12 % z). The summary vignette's claim that profileTLS and bayesTLS "agree across every row" (`case-study-summary.Rmd:246`) **overstates** this; the leaf-PSII vignette itself reports the gap honestly. Plausibly a prior-vs-prior-free effect on the Beta `phi`/asymptotes for a continuous endpoint.
- **Intervals are a different story from point estimates.** For **beta-binomial** data the profile CIs **under-cover**: cached coverage studies show CTmax/z coverage of ~0.84–0.88 at φ=50 and **collapsing to ~0.70 at φ=200**, because the profiled-out `phi` runs to the binomial limit and narrows the profile. Wald and bootstrap stay near nominal. profileTLS now warns and falls back to Wald — but the blanket "profile coverage tracks nominal 95 %" prose is true only for the **binomial** family.
- **Reproducibility caveat:** the benchmark cache stores only bayesTLS + two-stage summaries; **profileTLS estimates are not cached** — they are re-fit live in vignettes/tests. The static comparison tables are hand-written and could silently drift from the code. They currently match (I had them re-fit live), but nothing pins them.

**Net:** results are the same for well-identified count-data **point estimates**; they are **not** the same for the Beta/continuous case, and the **interval calibration differs** for beta-binomial. And on shrimp specifically, "same results" is comparing different (corrected vs corrupt) data.

---

## 4. How much time is saved?

- Cached timing (`timing_results.rds`, shrimp, one machine, 2026-06-18): **bayesTLS 5.01 s** (4 chains × 4000, *post-compile*), **two-stage 1.36 s**, **profileTLS ≈ 1 s** (profile path; Wald path near-instant; fits alone 6–87 ms in `performance-study.R`).
- **Claimed:** "comparable to two-stage and roughly an order of magnitude faster than the Bayesian MCMC fit."
- **Fairness caveats (so the speedup isn't overstated):**
  - bayesTLS was timed with **no `cores=` set** (`timing-study.R:41-44`) → likely 4 serial chains; with parallel chains bayesTLS could be ~4× faster, shrinking the gap.
  - The 5.01 s **excludes the one-time Stan compile** (favourable to bayesTLS, and disclosed).
  - profileTLS was timed live on a different machine/run than the cached comparators; n = 1, one seed (meta says "indicative, not a benchmark").
- **Honest reading:** the likelihood path is **genuinely much faster** — milliseconds-to-~1 s vs seconds-to-minutes (the latter including the unavoidable first-fit Stan compile). "Order of magnitude" holds for the **Wald** path; for the **profile** path against a serial 4-chain bayesTLS it's closer to **~5×**. Directionally correct, specific multiplier generous.

---

## 5. profileTLS errors & uncertainties (ranked)

**Errors / invalid steps**
- **E1 — Silent non-convergence (verified in `fit_engine.R`).** Non-convergence (`opt$convergence ≠ 0`) and non-positive-definite Hessian are *stored* (`convergence$code`, `pdHess`) but **no warning is emitted at fit time** (it only aborts if *both* optimisers error). A user who goes straight to `confint(method="wald")` or `$estimates` gets numbers from a bad fit with no signal. `confint(method="profile")` does react to `pdHess=FALSE`; the **Wald path does not.**
- **E2 — `up` confidence interval can fall outside [0,1].** `extract.R:223-230` builds the `up` interval as `up ± z·SE` on the natural scale (delta method) with no clamp; near the common `up≈1` saturation this routinely exceeds 1. Shipping an out-of-range CI for a bounded asymptote is a concrete defect (the bootstrap path stays in-range).
- **E3 — BFGS fallback adopts a non-converged optimum** (`fit_engine.R:64-74`) without flagging, compounding E1.
- **E4 — beta-binomial profile CIs under-cover** (see §3). Mitigated by a new weak-φ→Wald fallback, but the trigger `SE(φ)/φ̂ > 1` is an **uncalibrated magic threshold**, and "Wald stays calibrated when φ is weakly identified" is asserted more confidently than the in-repo evidence supports. It *does* warn (not silent).

**Fragilities / uncertain**
- **F1 — Boundary calibration not handled.** The "boundary" warning only fires when the deviance dips below 0 (reported MLE not the optimum), **not** the classic MLE-on-boundary case; `phi` at its boundary is still profiled with a plain `qchisq(·,1)` cutoff rather than the 50:50 χ²₀:χ²₁ mixture. (Variance components are kept on Wald, so moot for those.)
- **F2 — `uniroot` edge cases.** Multimodal profiles are *detected and warned* but the solve can still return the inner root (too-narrow interval); slow-but-closing profiles can be over-called "open." Both mostly rescued by the bootstrap fallback.
- **F3 — Random effects with few groups.** σ biased low (rel-bias −0.30 at 3 groups → −0.04 at 30) and CTmax Wald coverage 0.73→0.93; documented, with a "<~8 groups → prefer bayesTLS" advisory. RE *are* integrated out via Laplace correctly when profiling fixed targets.
- **F4 — Bootstrap is percentile-only** (no basic/BCa), so it inherits estimator skew; the construction-scale + back-transform choice keeps intervals in-range and exactly equivariant (good), and seed/parallel handling is correct and reproducible.
- **F5 — Brittle positional matching** of grouped coordinates to `name_map` (correct today, would silently mis-assign CIs if the C++ parameter order changed); several uncited magic constants (`3.5·se`, `floor_n=20`, `rel_se>1`).

**Minor:** a dead extra refit (`dev_at_hat`), a stale "bootstrap doesn't carry the random block" comment that contradicts the implemented RE bootstrap.

**What's genuinely good:** D = 2(ℓ̂−ℓₚ) and the χ²₁ cutoff are correct; equivariance `ci_z == exp(ci_log_z)` is achieved by construction (profile the unconstrained coordinate, transform endpoints); open sides return NA + warning (never a fabricated bound); nuisance parameters (incl. φ and RE) are genuinely re-optimised; the bootstrap DGP is consistent with the fitted family/φ/RE. The 12-warning identifiability story is a real value-add bayesTLS lacks.

---

## 6. Scope has grown beyond the package's own SPEC

`SPEC.md` §7 lists as **v0.1 non-goals**: Beta/continuous, time-to-event, multi-trait, heat-injury/repair, temp effects on low/up/k, random effects. The actual code now implements **Beta/continuous** (`family_code==2`, snowgum), **random effects** (CT/logz/low/logk intercepts), **heat injury + repair** (`heat_injury.R`), and shape covariates. This is ambitious and mostly well-built, but: (a) the SPEC/ROADMAP/known-limitations are likely **stale relative to the code**; and (b) the post-v0.1 additions are exactly where the weak spots concentrate (Beta snowgum mismatch, beta-binomial coverage, RE few-group bias). Treat the count-data binomial path as mature and the Beta/RE paths as newer.

---

## 7. Recommendations

**For us (bayesTLS) — highest priority:**
1. Fix the shrimp data bug (`make_datasets.R`: `deaths = round(prop*N)`), rebuild `shrimp_lethal.rda`, refit, re-render, and report whether the shrimp `z`/`CTmax`/`T_crit` change. This affects the manuscript.
2. Add a `data/` integrity check (the corrupt column collapses to {0,1} — a one-line assertion would have caught it).
3. Credit Shinichi/Nakagawa's team for surfacing it.

**For profileTLS (if asked to advise — not our package to edit):**
- Surface convergence/`pdHess` as a warning at fit time (E1/E3); clamp the `up` CI to [0,1] (E2).
- Default beta-binomial CIs to Wald/bootstrap or warn prominently given the profile under-coverage (E4); calibrate or document the `rel_se>1` threshold.
- Cache profileTLS benchmark estimates (prevent silent drift of the static tables).
- Soften two over-claims: "agree across every row" (snowgum doesn't) and "order of magnitude faster than MCMC" (≈5× for the profile path vs a serial-chain bayesTLS); time bayesTLS with `cores=4` for a fair comparison.
- Sync SPEC/ROADMAP/known-limitations to the actual (expanded) scope.

**Bottom line:** the statistics and calculations are sound and genuinely equivalent to bayesTLS for count-data point estimates, delivered much faster; the honest caveats are the Beta case, beta-binomial interval coverage, a few un-surfaced failure modes — and, most importantly for us, the shrimp data bug it exposed in our own package.
