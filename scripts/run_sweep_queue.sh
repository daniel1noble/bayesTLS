#!/bin/bash
# Sequential dispatcher for the three sensitivity sweeps:
#   A) u_beta1 grid (temperature slope on upper asymptote)
#   B) u_0 grid     (upper-asymptote level at T_bar)
#   C) design x n_reps grid (sparsity)
#
# All cells run at N = 500 sims, 5 PSOCK workers. Per-sim RDS files plus the
# aggregate (per_sim, meta, draws, summary, diffs) are saved per cell under
# the standard output/sim_twostage/ tree. Each cell logs to /tmp/sim_<label>.log;
# the queue itself logs every transition to /tmp/sim_sweeps.log.

set -e

WORKERS=5
# n_sims is set per scenario block below (1000 for Scen 1-5; 500 for sweeps).
# Each queue_cell call passes --n_sims explicitly via "$@".

# Pre-flight: build a DGP feasibility table for every planned cell. Each cell
# calls sim_twostage_truth() with the cell's parameters; the bounds check in
# compute_ols_truth() rejects DGPs whose u_T, ell_T, or k_T leave the
# (0, 0.5, 1) feasibility envelope at any design temperature. The check
# itself costs a few seconds total — much cheaper than discovering at hour 3
# that a sweep cell silently generates NaN survival probabilities.
preflight () {
  local label="$1"; shift
  # Pre-flight checks per cell:
  #   1. Truth function accepts the cell's CLI args (DGP feasibility guard
  #      in compute_ols_truth fires on infeasible u/ell/k).
  #   2. sim_twostage_dataset() runs with the cell's truth (catches NaN in
  #      y from edge-of-feasibility DGPs that the truth guard let through).
  #   3. Emit a hash of the resulting y vector to /tmp/sim_preflight_hashes
  #      so the queue can assert pairwise distinctness across cells.
  #      This is the check that would have caught the 2026-05-15 silent-
  #      baseline bug, where every cell scored against a shifted truth but
  #      the data stayed baseline. See feedback_sim_preflight.md.
  Rscript -e "
    source('scripts/sim_twostage_helpers.R')
    args <- commandArgs(trailingOnly = TRUE)
    parse_num <- function(key) {
      i <- match(paste0('--', key), args)
      if (is.na(i) || i == length(args)) NULL else as.numeric(args[i + 1])
    }
    parse_chr <- function(key, default) {
      i <- match(paste0('--', key), args)
      if (is.na(i) || i == length(args)) default else args[i + 1]
    }
    parse_int <- function(key, default) {
      i <- match(paste0('--', key), args)
      if (is.na(i) || i == length(args)) default else as.integer(args[i + 1])
    }
    truth <- sim_twostage_truth(
      dgp       = parse_chr('dgp',      'baseline'),
      u_0       = parse_num('u_0'),
      ell_0     = parse_num('ell_0'),
      u_beta1   = parse_num('u_beta1'),
      ell_beta1 = parse_num('ell_beta1'),
      k_beta1   = parse_num('k_beta1'),
      design    = parse_chr('design',   'full'),
      family    = parse_chr('family',   'beta_binomial')
    )
    # Resolve n_reps the same way sim_twostage_bias.R does.
    n_reps_cli <- parse_int('n_reps', NA_integer_)
    if (!is.na(n_reps_cli)) {
      n_reps <- n_reps_cli
    } else {
      scen <- parse_chr('scenario', 'n3')
      n_reps <- c(n3 = 3L, n5 = 5L)[[scen]]
    }
    # Generate one dataset with a fixed preflight seed (independent of the
    # master sim seed). If the CLI overrides are properly threaded, two
    # cells with different truths will produce different y vectors.
    d <- sim_twostage_dataset(n_reps = n_reps, seed = 1L, truth = truth)
    if (any(is.na(d\$y)))
      stop(sprintf('preflight: %s produces %d/%d NA y values — DGP infeasible',
                   '$label', sum(is.na(d\$y)), nrow(d)))
    h <- substr(rlang::hash(d\$y), 1, 12)
    cat(sprintf('  OK %s  z_true=%.4f  CTmax_1hr_true=%.4f  n_reps=%d  data_hash=%s\n',
                '$label', truth\$z_true, truth\$CTmax_1hr_true, n_reps, h))
    # Append the hash to the cross-cell hash file for queue-level distinctness check.
    cat(sprintf('%s\t%s\n', '$label', h),
        file = '/tmp/sim_preflight_hashes', append = TRUE)
  " "$@"
}

run_cell () {
  local label="$1"
  shift
  echo "[$(date '+%H:%M:%S')] === $label start ==="
  Rscript scripts/sim_twostage_bias.R \
      --label "$label" \
      --workers "$WORKERS" \
      "$@" \
      > "/tmp/sim_${label}.log" 2>&1
  echo "[$(date '+%H:%M:%S')] === $label done ==="
}

# ---- Collect every planned cell and pre-flight before running anything ----
# (This is a dispatcher pattern: build the queue first, validate, then run.)
declare -a CELL_LABELS CELL_ARGS
queue_cell () {
  CELL_LABELS+=("$1"); shift
  CELL_ARGS+=("$*")
}

# ============================================================================
# Scenarios are described in notes/2026-05-13-two-stage-bias-sim-methods.qmd.
# Naming convention: scenN_<name>_<rep label>. Scenarios 1-5 run at
# n_sims = 1000; sensitivity sweeps (6-8) at n_sims = 500.
# ============================================================================
NSIMS_MAIN=1000   # Scenarios 1-5
NSIMS_SWEEP=1000  # Scenarios 6-8 — standardised to match Scen 1-5 (was 500)

# -- Scenario 1: Strict equivalence baseline (binomial DGM) ------------------
# u=0.999, ell=0.001 (asymptotes at the 4PL bound), only midpoint varies with T.
# Both two-stage variants AND the joint 4PL are correctly specified.
for nr in 3 5; do
  queue_cell "scen1_strict_eq_n${nr}" --n_sims "$NSIMS_MAIN" \
    --dgp baseline --family binomial --u_0 0.999 --ell_0 0.001 --n_reps "$nr"
done

# -- Scenario 2: Likelihood misspec only (beta-binomial DGM, same shape) -----
# Same truth shape as Scen 1; only the DGM family changes (now BB with phi=5).
for nr in 3 5; do
  queue_cell "scen2_lik_misspec_n${nr}" --n_sims "$NSIMS_MAIN" \
    --dgp baseline --family beta_binomial --u_0 0.999 --ell_0 0.001 --n_reps "$nr"
done

# -- Scenario 3: Heat lowers max survival (asymmetric u drift) ---------------
# u drops 1pp/°C; ell constant at 0.05.
for nr in 3 5; do
  queue_cell "scen3_heat_lowers_u_n${nr}" --n_sims "$NSIMS_MAIN" \
    --dgp asym_u --family beta_binomial --n_reps "$nr"
done

# -- Scenario 4: Asymptotes compress (symmetric u + ell drift) ---------------
# u drops 1pp/°C, ell rises 1pp/°C.
for nr in 3 5; do
  queue_cell "scen4_compress_n${nr}" --n_sims "$NSIMS_MAIN" \
    --dgp sym_ul --family beta_binomial --n_reps "$nr"
done

# -- Scenario 5: Heat sharpens curve (k varies with T) -----------------------
# k goes from 7 at T=30 to 9 at T=38; asymptotes constant.
for nr in 3 5; do
  queue_cell "scen5_sharpen_n${nr}" --n_sims "$NSIMS_MAIN" \
    --dgp varying_k --family beta_binomial --n_reps "$nr"
done

# -- Scenario 6: Sensitivity — strength of u drift ---------------------------
# u_beta1 grid bounded by |u_beta1| < 0.020 (feasibility: u_T < 1 at T=30 with
# u_0 = 0.92). Deepest cell is -0.019 (~0.4pp buffer).
for ub in -0.005 -0.010 -0.015 -0.019; do
  ub_lab=$(awk -v v="$ub" 'BEGIN{ v=-v; printf "m%03d", v*1000 + 0.5 }')
  queue_cell "scen6_ub_${ub_lab}" --n_sims "$NSIMS_SWEEP" \
    --dgp baseline --family beta_binomial --u_beta1 "$ub" --n_reps 5
done

# -- Scenario 7: Sensitivity — level of upper asymptote ----------------------
for u0 in 0.99 0.95 0.85 0.75 0.65; do
  u0_lab=$(awk -v v="$u0" 'BEGIN{ printf "%03d", v*100 + 0.5 }')
  queue_cell "scen7_u0_${u0_lab}" --n_sims "$NSIMS_SWEEP" \
    --dgp baseline --family beta_binomial --u_0 "$u0" --n_reps 5
done

# -- Scenario 8: Sensitivity — design intensity ------------------------------
for cfg in "full 1" "full 3" "full 5" "sparse 1" "sparse 3" "sparse 5"; do
  read d n <<< "$cfg"
  queue_cell "scen8_${d}_n${n}" --n_sims "$NSIMS_SWEEP" \
    --dgp baseline --family beta_binomial --design "$d" --n_reps "$n"
done

# ---- Optional scenario filter (CLI args) -----------------------------------
# Run only cells whose label starts with one of the supplied prefixes.
# Examples:
#   ./scripts/run_sweep_queue.sh                # all 23 cells
#   ./scripts/run_sweep_queue.sh scen1 scen2 scen3   # only Scenarios 1-3 (6 cells)
if [ "$#" -gt 0 ]; then
  KEEP_LABELS=()
  KEEP_ARGS=()
  for i in "${!CELL_LABELS[@]}"; do
    for prefix in "$@"; do
      case "${CELL_LABELS[$i]}" in
        "$prefix"*)
          KEEP_LABELS+=("${CELL_LABELS[$i]}")
          KEEP_ARGS+=("${CELL_ARGS[$i]}")
          break
          ;;
      esac
    done
  done
  CELL_LABELS=("${KEEP_LABELS[@]}")
  CELL_ARGS=("${KEEP_ARGS[@]}")
  echo "[$(date '+%H:%M:%S')] Filtered queue to ${#CELL_LABELS[@]} cells matching prefixes: $*"
fi

# ---- Pre-flight pass: validate every planned cell BEFORE any long compute --
echo "[$(date '+%H:%M:%S')] Pre-flighting ${#CELL_LABELS[@]} cells ..."
rm -f /tmp/sim_preflight_hashes
for i in "${!CELL_LABELS[@]}"; do
  # shellcheck disable=SC2086 # intentional word-split of stored arg string
  preflight "${CELL_LABELS[$i]}" ${CELL_ARGS[$i]}
done

# Cross-cell distinctness check. Any two cells whose CLI args differ MUST
# generate different data on the same seed. A duplicate hash means a CLI
# flag isn't being threaded through to data generation (the 2026-05-15
# bug class — see feedback_sim_preflight.md).
#
# Known-equivalent pairs (allowed duplicates): scen3_heat_lowers_u_n5
# and scen6_ub_m010 share identical DGP arguments — Scen 6's u_beta1 =
# −0.010 cell IS the Scen 3 (asym_u) DGP at n_reps = 5. Listed below
# (semicolon-separated, order-tolerant) as an explicit whitelist; any
# other duplicate is treated as a real bug.
ALLOWED_DUPS="scen3_heat_lowers_u_n5,scen6_ub_m010;scen6_ub_m010,scen3_heat_lowers_u_n5"

if [ -f /tmp/sim_preflight_hashes ]; then
  # BSD-portable wide-uniq: group rows by the hash column (field 2),
  # collect labels per hash, and report any hash with >1 cell label.
  dup_pairs=$(awk '{ hashes[$2] = (hashes[$2] ? hashes[$2] "," $1 : $1) }
                   END { for (h in hashes) if (index(hashes[h], ",") > 0)
                           print hashes[h] "\t" h }' \
              /tmp/sim_preflight_hashes)
  if [ -n "$dup_pairs" ]; then
    # Filter out duplicates that match the whitelist (semicolon-separated).
    unexpected=$(echo "$dup_pairs" | awk -F'\t' -v allowed="$ALLOWED_DUPS" '
      BEGIN { n = split(allowed, a, ";"); for (i = 1; i <= n; i++) wl[a[i]] = 1 }
      { if (!(wl[$1])) print }')
    if [ -n "$unexpected" ]; then
      echo ""
      echo "  ERROR: unexpected duplicate data hashes across pre-flight cells."
      echo "  Cells producing identical data despite different CLI args:"
      echo "$unexpected" | awk -F'\t' '{print "    " $1 " (hash=" $2 ")"}'
      echo "  A CLI flag is not threaded through to data generation — see"
      echo "  feedback_sim_preflight.md."
      exit 1
    fi
    # Whitelisted duplicates are fine — log them for transparency.
    echo "$dup_pairs" | awk -F'\t' '{print "  NOTE: expected duplicate "
                                          $1 " (hash=" $2 ") — by-design overlap."}'
  fi
fi
echo "[$(date '+%H:%M:%S')] Pre-flight passed (cross-cell distinctness check OK)."

# ---- Dispatch ------------------------------------------------------------
N_CELLS=${#CELL_LABELS[@]}
QUEUE_START=$(date +%s)
for i in "${!CELL_LABELS[@]}"; do
  CELL_NUM=$((i + 1))
  echo ""
  echo "[$(date '+%H:%M:%S')] === Queue progress: cell ${CELL_NUM} / ${N_CELLS} (${CELL_LABELS[$i]}) ==="
  # shellcheck disable=SC2086 # intentional word-split of stored arg string
  run_cell "${CELL_LABELS[$i]}" ${CELL_ARGS[$i]}
  # Per-cell summary: the driver already prints mcse_summary + diag_summary
  # at the end of its log; tail just those headlined sections so the queue
  # log gives a heartbeat without dragging in the full per-sim chatter.
  echo "[$(date '+%H:%M:%S')] --- Cell ${CELL_LABELS[$i]} summary (last 80 lines of cell log) ---"
  tail -n 80 "/tmp/sim_${CELL_LABELS[$i]}.log"
  ELAPSED=$(( $(date +%s) - QUEUE_START ))
  H=$((ELAPSED / 3600)); M=$(((ELAPSED % 3600) / 60))
  if [ "$CELL_NUM" -lt "$N_CELLS" ]; then
    AVG=$((ELAPSED / CELL_NUM))
    REMAIN=$(( AVG * (N_CELLS - CELL_NUM) ))
    RH=$((REMAIN / 3600)); RM=$(((REMAIN % 3600) / 60))
    echo "[$(date '+%H:%M:%S')] === ${CELL_NUM}/${N_CELLS} cells complete | elapsed ${H}h${M}m | est remaining ${RH}h${RM}m ==="
  else
    echo "[$(date '+%H:%M:%S')] === All ${N_CELLS} cells complete | total elapsed ${H}h${M}m ==="
  fi
done

echo "[$(date '+%H:%M:%S')] All sweeps complete"
