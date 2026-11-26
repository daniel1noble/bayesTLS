# Provenance of the two-stage-bias simulation results

The simulation compares the classical two-stage pipeline with the joint 4PL
model across 29 scenario-by-design cells of 1,000 simulated datasets each. Its
raw per-simulation files (`output/sim_twostage/raw/`, on OSF as the `full`
tier) and the summaries the manuscript and supplement read (`per_sim_*`,
`summary_*`, `diffs_*`, `diag_*`, `meta_*`; the OSF `results` tier) come from
two runs.

## The June 2026 run (all 29 cells)

- `scripts/simulations/run_simulations.R`, 11-13 June 2026, on R 4.5 (macOS,
  arm64) with CmdStan 2.38.0 and the bayesTLS code of mid-June 2026, whose
  `extract_tdt()` read z, CTmax_1hr and T_crit from the joint fits.
- Package versions were not recorded with the results. The R 4.5 library on
  the machine that ran it holds glmmTMB 1.1.14, TMB 1.9.17, Matrix 1.7-4,
  brms 2.23.0 and cmdstanr 0.9.0.

## The November 2026 re-run (14 cells, joint 4PL rows only)

- Why: `extract_tdt()` inverted each draw's absolute-threshold LT curve on a
  fixed grid (`min(T) - 2` to `max(T) + 2` C) and silently dropped draws whose
  CTmax_1hr crossing lay outside it, so the "Joint 4PL (absolute)" CTmax_1hr
  summaries of those cells came from truncated posteriors.
- Code: `scripts/simulations/rerun_abs_truncation.R`, with `extract_joint()`
  in `sim_functions.R` replacing `extract_tdt()` (see its documentation for
  exactly what changed). `check_rerun_scope.R` first checked, on the first 6
  simulations of every cell, that the re-run reproduces the June joint fits:
  the relative arm to 4.9e-10 and the absolute z to 8.2e-9.
- Software: R 4.6.2; CmdStan 2.38.0 (2.40.0 does not reproduce the June
  draws); a frozen bayesTLS (package code of commit a4169a4) in a private
  library with glmmTMB 1.1.15.2, TMB 1.9.25 and Matrix 1.7-6; brms 2.23.0,
  cmdstanr 0.9.0 and posterior 1.7.0.
- Merge rule (`merge_rerun()`): in a re-run simulation, the joint 4PL rows
  (relative and absolute) come from the re-run; the two-stage rows, the
  two-stage success flags and the sampler diagnostics come from the June
  file. The two-stage rows are kept because with glmmTMB 1.1.15.2 they differ
  from June's in 30 of the 174 scope-check simulations (by more than 1e-6; 18
  by more than 1e-3, at most 13.9). Cells that were not re-run are the June
  files unchanged.

## Reproducing

- Joint fits: CmdStan 2.38.0 (`run_simulations.R` pins it).
- Two-stage rows: the June glmmTMB version (1.1.14 on R 4.5, above);
  other versions change some beta-binomial fits.
- Re-run cells: `rerun_abs_truncation.R` documents its library, CmdStan and
  folders, and keeps the June raw files in `raw_pre_rerun_2026-11-26/`.
