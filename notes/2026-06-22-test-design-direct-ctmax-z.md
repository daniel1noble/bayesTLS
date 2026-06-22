# Test plan: direct-CTmax/z redesign — audit of existing tests + new test design

**Status:** test design / for review before implementation. **Date:** 2026-06-22.
Synthesised from a 14-agent audit+design+adversarial-critique pass. Builds on the
design note (`2026-06-22-direct-ctmax-z-package-design.md`) and the prototype
validation (`2026-06-22-direct-ctmax-z-prototype-validation.qmd`).

## 0. The headline finding (why the obvious test plan is not enough)

The existing suite has a deliberate **fast/skip split**: formula-string, prior-object
and planted-draws tests run on every `devtools::test()`/CI run; every test that
*fits* a model is gated behind `skip_unless_brms()` (runs only when
`RUN_BRMS_TESTS=true` **and** cmdstanr is installed). Default CI never compiles Stan.

The naive test plan put **all** direct-mode correctness checks (silent-NA guard,
threshold conversion, equivalence, divergence contrast) into the gated tier. So on
every normal push **zero tests would verify the redesign is correct** — only that
the formula strings parse. The single most dangerous regression we found — a direct
fit makes `extract_tdt()` return **all NA, silently** (verified: the relative CTmax
branch computes `Tbar + (target - b_mid_Intercept)/b_mid_temp_c`; on a direct fit the
missing names read as 0, so `Tbar + x/0 = Inf` for every draw → `filter(is.finite)` →
0 rows → NA) — would have **no default-CI guard at all**.

So the plan below is reorganised around one principle: **the correctness-critical
checks must be fast, deterministic and run by default.** That is achievable, but it
forces three implementation decisions (below).

## 1. Decisions the tests force

> **LOCKED 2026-06-22 (Daniel approved all five):** D-T1 rewire the three helpers to
> read coefficients (CTmax/z are estimated coefficients; `posterior_linpred` only for
> continuous covariates on CTmax/z) · D-T2 add a committed `seed` arg to
> `tls()`/`extract_tdt()` · D-T3 factor `resolve_shape()` as a named internal function
> · D-T4 encode the `lower/upper`→`bounds` policy as a test · D-T5 one gitignored +
> auto-rebuilt grouped direct fixture for gated tests (migrate `sim_fit_beta_small.rds`
> to match), case-study golden numbers go to a scripted (CI-rendered) validation gate,
> not testthat.

- **D-T1 — rewire the three helpers to READ COEFFICIENTS, not delegate to
  `posterior_linpred`/`tls()`.** The design note offered "route through
  `posterior_linpred` or delegate to `tls()`." The audit shows that choice makes the
  direct-mode tests (a) require a real Stan fixture (a fake `as_draws_df` cannot
  answer `posterior_linpred`), and (b) tautological (`extract_tdt == tls` shares the
  same `C(T)` code). **Recommendation:** for the common intercept / factor-level case,
  `extract_tdt()`/`predict_heat_injury()`/`tdt_parameter_table()` read
  `z = exp(b_logz_<coef>)`, `CTmax = temp_mean + b_CTmaxdev_<coef>` off the
  `as_draws_df` (the prototype §9 "CTmax/z ARE coefficients"); keep `posterior_linpred`
  only as the general path for *continuous covariates* on CTmax/z. This makes the
  silent-NA / dose / C(T) guards fast, deterministic, brms-free **and** non-tautological
  (truth-anchored). (`tls()` stays as-is — its grid-slope path already works.)
- **D-T2 — add a committed `seed` argument to `tls()`/`extract_tdt()`**, threaded to the
  `runif()` site(s) used for the rate-multiplier T_crit (`R/tls.R:168`,
  `R/extract_tdt.R:644`). These are currently **unseeded** → T_crit is non-reproducible
  run-to-run, in *midpoint mode too* (a standing CLAUDE.md §4 RNG-reproducibility gap).
  Without this, any T_crit equality test is flaky; with it, determinism is testable.
- **D-T3 — factor `resolve_shape()` as a NAMED internal function** (not an inlined
  closure), so the inheritance truth-table and its `gsub` edge cases (RE terms with
  spaces `(1 | Date)`, multi-term cell-means, treatment vs cell-means) are unit-testable
  string-in/string-out.
- **D-T4 — `lower`/`upper` → `bounds` migration policy** must be *encoded*, not just
  decided: either `expect_error(fit_4pl(std, lower = 0.85), "lower.*(deprecated|bounds)")`
  (if removed — guards the silent-swallow-by-`...` where a sublethal model misfits at
  bounds (0,1)), or `expect_identical(make_4pl_formula(bounds = c(0.85,1)),
  make_4pl_formula(lower = 0.85, upper = 1))` (if aliased).
- **D-T5 — fixture strategy:** ONE consolidated gated direct fixture (see §4);
  case-study golden numbers move to a *scripted* gate, not testthat (see §5); resolve the
  committed-vs-gitignored fixture inconsistency.

## 2. Existing-test audit (what breaks, what guards backward-compat)

Verdicts (full table in the workflow output). The redesign is **additive**, so almost
nothing *breaks* — most existing tests cover midpoint mode (unchanged) and thereby
become the **backward-compatibility guard**:

| verdict | files | action |
|---|---|---|
| backward-compat-guard | `test-fit_4pl.R` | KEEP all midpoint structural assertions (mid nlpar present, `temp_effects` toggle, RE-on-mid, "mid must always carry" error). Update only the `lower/upper`→`bounds` calls per D-T4. |
| needs-update | `test-priors.R`, `test-diagnostics.R` | extend for direct-mode priors / `tdt_parameter_table` direct path; keep midpoint assertions. |
| should-be-extended | `helper-skip.R`, `helper-simulate.R` | add `load_fixture_workflow_direct()` + a grouped/asymmetric/optionally-minutes simulator and a per-mode (relative AND absolute, per-group) truth oracle. |
| safe (≈18 files) | extract_tdt, tls, predict_*, accessors, derive_z, etc. | unchanged; they pin the midpoint path. **Keep the strict `1e-8` identity in `test-extract_tdt.R:73` — do not loosen it to a "byte-identical" check.** |

**Infrastructure (must match for new tests):** fits go through
`load_fixture_workflow*()` → `fit_4pl(file = fixtures/sim_fit_small.rds)` (brms file-cache),
all gated by `skip_unless_brms()`; both committed fixtures are **midpoint**, so they
cannot exercise direct mode. Fast tests build fake `as_draws_df` of known coefficients
(`fake_workflow`, `fake_hi_workflow`) or call internal helpers on synthetic draws — all
seeded/deterministic.

## 3. New tests — Tier 1: FAST, deterministic, run on every CI push

These are the bulk and the real protection. No Stan.

**Formula resolution / inheritance (string-level; require D-T3 so `resolve_shape` is testable):**
- Supplying `ctmax`/`z` triggers direct mode: estimated nlpars become `CTmaxdev`+`logz`, `mid` becomes `nlf`-derived; supplying **neither** reproduces today's midpoint formula **exactly**.
- Inheritance truth-table on `resolve_shape()`: intercept-only → `~ temp_c`; cell-means `0+G` → `0+G+temp_c:G`; treatment `G` → `G*temp_c`; RE **not** inherited; explicit wins; inheritance source `ctmax` (fallback `z`); pooled-z (`z` omitted → `logz ~ 1`). *Assert on the brms-parsed coefficient set (`get_prior()`/`brmsterms()`), not the deparsed string* — so term reordering (`temp_c:G` vs `G:temp_c`) is tolerated but the actual coefficients are pinned, and the formula round-trips through brms without error.
- `threshold = "absolute"` adds the `-C(T) = (1/k)log((up-p)/(p-low))` term to the `mid` line; relative omits it; non-default `p` flows in; `bounds = c()` bakes the same asymptote constants as the old `lower/upper`.

**Priors (object-level):**
- **Centred lowraw/upraw/logk Intercept priors are byte-identical between midpoint and direct mode** (`expect_identical`) — *this is the only reproducible guard for the 191-vs-0 divergence finding* (see §5 for why the sampling test cannot be).
- **Frozen-snapshot** midpoint priors: assert the *literal* expected strings to 6 dp (`lowraw` Intercept `normal(-3.227262, 1)`, `logk` `normal(0.693147, 1)`, `mid` `normal(<median logd>, 1.5)`, sds, classes, nlpars) — hard-coded, **not** re-derived from `make_4pl_priors` internals (re-deriving is a tautology; only frozen literals catch a shared refactor that flattens *both* modes).
- Direct priors: `CTmaxdev`/`logz` b-priors present, **no `mid` prior**; cell-means `ctmax` → one `CTmaxdev` prior per level; RE on ctmax/z → sd prior on `CTmaxdev`/`logz` **only** (not on shape); `get_prior()` completeness (every estimated coef has a prior); phi present for beta/beta-binomial, absent for binomial.

**Planted-draws / fake-`as_draws_df` (require D-T1 coef-read path) — the correctness guards:**
- **Silent-NA guard (highest value):** build a fake direct `as_draws_df` (`b_lowraw_Intercept`, `b_upraw_Intercept`, `b_logk_Intercept`, `b_CTmaxdev_Intercept`, `b_logz_Intercept`, *no* `b_mid_*`), call `extract_tdt()`, assert **both** `z$summary$z_median` **and** `CTmax$summary$temp_median` are **finite** (the prototype probe only checked z; the CTmax=Inf→NA path is the one I verified), with `z == median(exp(b_logz))` and `CTmax == temp_mean + median(b_CTmaxdev)`.
- **C(T) arithmetic, independent oracle:** plant known `low/up/logk/CTmaxdev/logz`, hand-compute `C(T)=(1/k)log((up-p)/(p-low))`, assert the rewired absolute conversion reproduces it. Negative case: sublethal `up<p` → `C(T)` is `NaN` (demonstrates the O2 guard is load-bearing). *Independent of `tls()` — catches a shared-sign bug that `extract_tdt==tls` cannot.*
- **`predict_heat_injury` dose oracle:** fake single-group direct draws, run on a planted trace, assert accumulated dose == `planted_dose_from_trace(z = exp(b_logz), CTmax = temp_mean + b_CTmaxdev)` (the existing analytic oracle, computed independently of how the function reads coefficients) — catches the `mid_int not found` crash **and** a wrong-sign/scale read.
- **`/exp(logz)` finiteness:** feed `logz ∈ [-5, 5]` into the `mid` expression, assert no `Inf/NaN` propagates (the funnel-blowup risk; a "healthy fixture" check passes vacuously).
- **`log10_tref` units:** assert the helper returns `log10(60)` for `t_ref = 1h, duration_unit = "minutes"` and `0` for `"hours"` — the only fast guard against the `z·1.78 °C` CTmax shift from a t_ref-as-model-units bug.

**Backward-compat / API (string/signature-level):**
- Frozen-snapshot midpoint `make_4pl_formula()` deparse; midpoint priors keep nlpar set `{lowraw,upraw,logk,mid}` and gain **no** `CTmaxdev`/`logz`.
- Exported signatures: new args are additive (NULL-defaulted), no existing arg dropped/required.
- Legacy-arg policy test per **D-T4**.

## 4. New tests — Tier 2: GATED brms, ONE consolidated direct fixture

`load_fixture_workflow_direct()` — a **single** small **grouped** direct relative fit
(2-group sim, `ctmax/z ~ 0 + grp`, no RE, gitignored + auto-rebuild like
`sim_fit_small.rds`), gated by `skip_unless_brms()`. It serves all of:

- **Truth-recovery (THE non-tautological correctness anchor):** assert the rewired
  `extract_tdt()` recovers simulator-truth `z` and `CTmax_1hr` within tolerance / inside
  the CrI (reference derived independently of the code path, à la `test-extract_tdt.R:6`).
  Catches log10_tref/units, logz-scaling and silent-NA bugs simultaneously.
- **O1 explicit:** `expect_no_error` + finite on `posterior_linpred(fit, nlpar = "mid")` —
  the direct check that the `nlf`-derived `mid` is reachable (don't rely on it transitively
  through `tls`).
- **Per-group coef reads:** `z == exp(b_logz_<level>)` for each level *by name* (use a large
  known z gap so a swapped-group mapping can't pass by coincidence).
- **Threshold converter, correct modes + difference:** on this relative fit,
  `extract_tdt(target_surv="relative") == coef-read == tls(mode="relative")`;
  `extract_tdt(target_surv="absolute") == tls(mode="absolute")`; **and assert the two
  DIFFER** (by analytic `C(T)`) — guards a no-op or transposed-mode converter. (Corrects the
  design note's repeated mis-statement that the relative coef-read equals `tls(mode="absolute")`.)
- **Unchanged-helper canary:** `tls`/`predict_survival_curves`/`derive_tdt_landscape`/`diagnose_tdt_fit` run on the direct fit and give sane output.
- **Accessor round-trip** on the direct `extract_tdt()` result.
- **Backward-compat coefficient-absence:** load the MIDPOINT fixture after the rewire, assert
  `extract_tdt` z still `== median(-1/b_mid_temp_c)` at **1e-8**, AND the fit has **no**
  `b_logz`/`b_CTmaxdev` and **does** have `b_mid_temp_c` (the negative assertion catches a
  rewire that leaks the direct branch into the default path).

T_crit on the fixture: assert only runif-invariant quantities or pass the **D-T2** seed;
never a bare point-T_crit equality.

## 5. New tests — Tier 3: scripted pre-merge gate (NOT testthat)

These cannot live in the fixture regime (case-study magnitudes need the heavy
`output/models/*.rds`; the funnel needs case-study geometry). Put them in a CI-rendered
`notes/*.qmd` gate, scripted (not "manual", which rots):

- **Case-study equivalence:** zebrafish direct z (2.09/1.76/2.28) & CTmax (39.77/41.34/39.76)
  match the midpoint reference to ~0.05 °C; Drosophila z (F 3.05/M 3.20) & **CI width** match.
- **Divergence contrast (191 vs 0):** on a fixture *known to exhibit the funnel*
  (small grouped sim with sparse hot cells), with documented expected counts — **not** asserted
  at tiny single-group scale where it gives 0/0 (false pass).
- **Grouped ABSOLUTE per-group truth (the inheritance-bias catch):** fit a 2-group absolute
  model with deliberately different per-group asymmetry twice — `up/low/k` **inherited**
  (→ `temp_c:G`) vs forced **shared** `+temp_c`. Assert the inherited fit's per-group absolute
  CTmax matches per-group truth **and differs** from the shared-slope fit. *This is the only
  test of rule-2's rationale; all relative golden tests are provably insensitive to it.*

## 6. Tests we DROP (tautological — they pass even if the bug is present)

- `coef-read == posterior_linpred(nlpar="logz")` — re-asserts brms equals its own coefficients; exercises none of the rewired code.
- `tls(mode="relative") == exp(b_logz)` — holds *by the formula's algebra* (`mid` slope is exactly `-1/exp(logz)`); tests `make_4pl_formula`, not the extractor; passes even if `extract_tdt` returns all-NA.
- `z == median(exp(b_logz))` as the *only* check — `f(x)==x` when the code reads `exp(b_logz)`. Keep it **only** alongside the independent truth-recovery anchor.
- "byte-identical before/after the rewire" framed as an in-checkout diff — there is no "before"; replace with the **frozen literal snapshots** (§3).
- The wide `T_crit ∈ [CTmax-3z, CTmax-z]` window — too loose to distinguish correct from mis-paired draws; a cheap `is.finite` + truth-in-CrI is better.

## 7. Resolved

1. **D-T1–D-T5**: all approved (see §1).
2. **grouped `predict_heat_injury`**: still to settle during implementation — either add a group selector (full feature) or scope to single-group with an explicit `expect_error` on grouped input. Decide when that function is rewired; flag to Daniel at that point.

## 8. Implementation status (2026-06-22)

- **Tier 1 (fast)**: `test-direct-formula.R`, `test-direct-priors.R`,
  `test-direct-extract.R` — formula/prior structure, C(T) oracle, silent-NA
  guards, all four absolute-fit threshold cells, unit aliases, seed determinism,
  print parameterisation. Run every CI pass.
- **Tier 2 (gated)**: consolidated fixture `load_fixture_workflow_direct()` is a
  **single-condition** relative fit (not the 2-group originally sketched in §4):
  `extract_tdt` is single-condition by design, so per-group reads are validated
  via `tls()` in the Tier-3 gate instead. `test-direct-fixture.R` covers
  truth-recovery, direct↔midpoint equivalence, the integration silent-NA guard,
  grouped-redirect for every single-condition helper, and seed reproducibility.
- **Tier 3 (scripted gate)**: `notes/2026-06-22-direct-ctmax-z-validation-gate.qmd`
  — grouped ABSOLUTE per-group truth (inheritance-bias, hard `stopifnot`) and the
  divergence contrast (centred vs flat, qualitative). **Case-study golden numbers
  are deliberately deferred** while the case studies are being restructured
  (`ms/case_studies_new.qmd`); re-point a case-study section into the gate once
  the new studies settle. Equivalence is covered meanwhile by the Tier-2 fixture
  and `notes/2026-06-22-direct-ctmax-z-prototype-validation.qmd`.
- **Adversarial review (2026-06-22)** of the rewire confirmed and fixed: the
  `tdt_loglt` C0/C_p conflation, `derive_z` / `derive_temperature_for_duration`
  direct wiring, the treatment-coded grouped guard (`tdt_is_grouped`), and the
  `fit_4pl` unit-alias desync — all with regression tests.
