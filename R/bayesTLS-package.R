#' bayesTLS: Joint Bayesian 4PL Models for Thermal Load Sensitivity
#'
#' Functions for fitting a single joint Bayesian four-parameter logistic
#' model to thermal-tolerance proportion data and extracting all classical
#' TLS quantities (z, CTmax at 1 hour, T_crit) and downstream predictions
#' (heat-injury accumulation, predicted survival under fluctuating
#' temperature traces, optional Sharpe-Schoolfield repair) from a single
#' fitted model.
#'
#' @keywords internal
"_PACKAGE"

# Suppress R CMD check NOTEs about non-standard evaluation (NSE) used by
# dplyr verbs, ggplot2 aes(), and brms formula coefficients. None of these
# are global variables in the package sense; they are column names referred
# to inside data-aware expressions.
utils::globalVariables(c(
  # standardized column names produced by standardize_data() and used by the
  # prediction / extraction pipeline:
  ".draw", ".epred",
  "assay_temp", "duration", "duration_model", "duration_out",
  "duration_lower", "duration_median", "duration_upper",
  "hi", "hi_lower", "hi_median", "hi_upper",
  "log10_duration", "log10_rate", "log10_t", "logd",
  "mortality", "mort_lower", "mort_median", "mort_upper",
  "n", "n_total", "n_surv",
  "p_true", "Parameter",
  "r_squared",
  "scenario", "slope_T",
  "survival", "survival_lower", "survival_median", "survival_upper",
  "survival_mean", "survival_se", "n_units",
  "surv_lower", "surv_median", "surv_upper",
  "T_at",
  "target_surv", "temp", "temp_c",
  "temp_lower", "temp_median", "temp_upper",
  "time_h", "value", "model", "quantity",
  "z", "z_lower", "z_median", "z_upper",
  # raw posterior columns referenced in helpers / extractors:
  "b_lowraw_Intercept", "b_upraw_Intercept",
  "b_logk_Intercept", "b_mid_Intercept", "b_mid_temp_c",
  # 4PL natural-scale parameters reconstructed inside dplyr::mutate chains:
  "low", "up", "k", "mid", "mid_int", "mid_temp",
  "ell", "u", "ellraw", "uraw", "logk",
  # NSE references used in plotting helpers:
  "T_c", "duration_mid", "y_ci",
  # additional NSE references in derive_tdt_parameters / plot_repair_rate:
  "CTmax", "repair_rate"
))
