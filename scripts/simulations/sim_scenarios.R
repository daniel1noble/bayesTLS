# ============================================================================
# Two-stage-bias simulation -- the scenario table (single source of truth).
# Sourced by run_simulations.R and rerun_abs_truncation.R. The row order sets
# each scenario's `index`, and so its seeds (MASTER_SEED + 1000 * index +
# sim_id): never reorder or delete rows; append new ones at the end.
# ============================================================================
# Each row is plain data: a label plus the sim_truth() arguments for that cell.
# NA in u_0 / ell_0 / u_beta1 means "keep the DGP preset". Edit a number to
# change a sweep value; add a row to add a cell.
#
#   scen1  strict-equivalence baseline (binomial DGM, asymptotes at the bound)
#   scen2  likelihood misspecification only (beta-binomial, same shape)
#   scen3  heat lowers max survival      (asymmetric u drift)
#   scen4  asymptotes compress           (symmetric u + ell drift)
#   scen5  heat sharpens the curve       (k varies with T)
#   scen6  sweep: strength of u drift    (u_beta1)
#   scen7  sweep: level of upper asymptote (u_0)
#   scen8  sweep: design intensity       (full/sparse x replication)
#   scen9  sweep: exposure-time window   (max duration 60/120/240/405 min) -- the
#          design question: how well does a shortened assay recover z & CTmax_1hr?
scenarios <- tibble::tribble(
  ~label,                   ~dgp,        ~family,         ~n_reps, ~design,   ~u_0,  ~ell_0, ~u_beta1,
  "scen1_strict_eq_n3",     "baseline",  "binomial",      3,  "full",    0.999, 0.001,  NA,
  "scen1_strict_eq_n5",     "baseline",  "binomial",      5,  "full",    0.999, 0.001,  NA,
  "scen2_lik_misspec_n3",   "baseline",  "beta_binomial", 3,  "full",    0.999, 0.001,  NA,
  "scen2_lik_misspec_n5",   "baseline",  "beta_binomial", 5,  "full",    0.999, 0.001,  NA,
  "scen3_heat_lowers_u_n3", "asym_u",    "beta_binomial", 3,  "full",    NA,    NA,     NA,
  "scen3_heat_lowers_u_n5", "asym_u",    "beta_binomial", 5,  "full",    NA,    NA,     NA,
  "scen4_compress_n3",      "sym_ul",    "beta_binomial", 3,  "full",    NA,    NA,     NA,
  "scen4_compress_n5",      "sym_ul",    "beta_binomial", 5,  "full",    NA,    NA,     NA,
  "scen5_sharpen_n3",       "varying_k", "beta_binomial", 3,  "full",    NA,    NA,     NA,
  "scen5_sharpen_n5",       "varying_k", "beta_binomial", 5,  "full",    NA,    NA,     NA,
  "scen6_ub_m005",          "baseline",  "beta_binomial", 5,  "full",    NA,    NA,     -0.005,
  "scen6_ub_m010",          "baseline",  "beta_binomial", 5,  "full",    NA,    NA,     -0.010,
  "scen6_ub_m015",          "baseline",  "beta_binomial", 5,  "full",    NA,    NA,     -0.015,
  "scen6_ub_m019",          "baseline",  "beta_binomial", 5,  "full",    NA,    NA,     -0.019,
  "scen7_u0_099",           "baseline",  "beta_binomial", 5,  "full",    0.99,  NA,     NA,
  "scen7_u0_095",           "baseline",  "beta_binomial", 5,  "full",    0.95,  NA,     NA,
  "scen7_u0_085",           "baseline",  "beta_binomial", 5,  "full",    0.85,  NA,     NA,
  "scen7_u0_075",           "baseline",  "beta_binomial", 5,  "full",    0.75,  NA,     NA,
  "scen7_u0_065",           "baseline",  "beta_binomial", 5,  "full",    0.65,  NA,     NA,
  "scen8_full_n1",          "baseline",  "beta_binomial", 1,  "full",    NA,    NA,     NA,
  "scen8_full_n3",          "baseline",  "beta_binomial", 3,  "full",    NA,    NA,     NA,
  "scen8_full_n5",          "baseline",  "beta_binomial", 5,  "full",    NA,    NA,     NA,
  "scen8_sparse_n1",        "baseline",  "beta_binomial", 1,  "sparse",  NA,    NA,     NA,
  "scen8_sparse_n3",        "baseline",  "beta_binomial", 3,  "sparse",  NA,    NA,     NA,
  "scen8_sparse_n5",        "baseline",  "beta_binomial", 5,  "sparse",  NA,    NA,     NA,
  "scen9_tmax_060",         "baseline",  "beta_binomial", 5,  "tmax060", NA,    NA,     NA,
  "scen9_tmax_120",         "baseline",  "beta_binomial", 5,  "tmax120", NA,    NA,     NA,
  "scen9_tmax_240",         "baseline",  "beta_binomial", 5,  "tmax240", NA,    NA,     NA,
  "scen9_tmax_405",         "baseline",  "beta_binomial", 5,  "tmax405", NA,    NA,     NA,
)
scenarios$index <- seq_len(nrow(scenarios))   # disjoint seed offset per scenario

