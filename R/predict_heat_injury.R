# Heat injury and predicted survival under a fluctuating temperature trace.
# Carries the model's full posterior through to HI(t) and S(t) trajectories
# with credible intervals. The algorithm:
#
#   1. Extract per-draw 4PL parameters (low, up, k, mid_int, mid_temp).
#   2. At each time step, compute LT_target(T) analytically from the 4PL,
#      giving damage rate = 1 / LT_target.
#   3. Optionally add a temperature-dependent repair rate (Sharpe-Schoolfield).
#   4. Forward-Euler integrate (damage - repair) across the trace, using each
#      interval's own width (dose starts at zero at the first time point).
#   5. Map cumulative dose back to predicted survival via the 4PL.

#' Per-draw 4PL parameters from a fitted workflow
#'
#' Extracts `(low, up, k, mid_int, mid_temp)` per posterior draw, suitable for
#' the analytical heat-injury machinery. Temperature slopes on `lowraw`,
#' `upraw`, and `logk` are intentionally ignored — the heat-injury integral
#' is evaluated under the classical assumption that the asymptotes and slope
#' are constant in T (only the midpoint shifts). If those temp slopes are
#' shrunk near zero by the data, this introduces no bias; otherwise it is the
#' same approximation the classical HI framework makes.
#'
#' @param workflow Fitted `bayes_tls`.
#' @param by Optional moderator column(s) for per-group parameters. `NULL`
#'   (default) uses the fit's moderators; a single-condition fit returns the
#'   ungrouped tibble. A grouped fit prepends the moderator column(s).
#' @return A tibble with `(.draw, low, up, k, mid_int, mid_temp)` columns
#'         (plus the moderator column(s) for a grouped fit), filtered to draws
#'         producing valid parameter values.
#' @examples
#' \dontrun{
#' pars <- extract_4pl_pars(wf)
#' head(pars)
#' }
#' @export
extract_4pl_pars <- function(workflow, by = NULL) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  # Constant-in-T 4PL parameters for the heat-injury integral, read by evaluating
  # posterior_linpred(nlpar=) at temp_c = 0 (low/up/k + midpoint intercept) and
  # temp_c = 1 (for the linear midpoint slope) — parameterisation- and
  # coding-agnostic, no coefficient-name parsing. The classical constant-shape
  # assumption is preserved: low/up/k are read at temp_c = 0 only; mid is the one
  # T-varying quantity (intercept + slope). For a direct absolute fit the nlf mid
  # already folds in the (constant-shape) asymmetry correction, so it stays linear.
  fit <- get_brmsfit(workflow)
  by  <- tdt_resolve_by(workflow, by)
  nd  <- tls_build_grid(fit$data, by = by, temp = "temp_c", temp_grid = c(0, 1))
  sp  <- tls_eval_subpars(fit, nd, workflow$meta$bounds, mode = "relative")

  per_group <- lapply(unique(nd$.grp), function(g) {
    gi <- which(nd$.grp == g)
    c0 <- gi[which(nd$temp_c[gi] == 0)]; c1 <- gi[which(nd$temp_c[gi] == 1)]
    out <- tibble::tibble(
      .draw    = seq_len(nrow(sp$mid)),
      low      = sp$low[, c0], up = sp$up[, c0], k = sp$k[, c0],
      mid_int  = sp$mid[, c0], mid_temp = sp$mid[, c1] - sp$mid[, c0]
    )
    out <- dplyr::filter(out,
                         is.finite(low) & is.finite(up) & is.finite(k) &
                         is.finite(mid_int) & is.finite(mid_temp) & k > 0 & up > low)
    if (!is.null(by)) out <- cbind(nd[gi[1], by, drop = FALSE], out, row.names = NULL)
    out
  })
  dplyr::bind_rows(per_group)
}

#' Analytical inverse 4PL: duration to reach a target survival at a given temperature
#'
#' For each `temp`, computes the exposure duration at which the 4PL gives
#' the chosen survival threshold. In `"relative"` mode the threshold is the
#' per-draw midpoint between the asymptotes, which collapses to the bare 4PL
#' `mid` parameter on the natural time axis. In absolute mode it is the
#' literal probability `survival_target`.
#'
#' @param temp           Numeric vector of temperatures (°C).
#' @param survival_target Either the literal probability to invert at (must
#'                       lie strictly between `low` and `up`), or the string
#'                       `"relative"` to use `(low + up)/2`.
#' @param low,up,k       Scalar (or per-temp) 4PL parameters.
#' @param mid_int,mid_temp Midpoint sub-model coefficients.
#' @param temp_mean      Centring temperature used by the model.
#' @return Numeric vector of durations in the model's time units. `NA` where
#'         the threshold is outside `(low, up)` in absolute mode.
#' @keywords internal
time_to_surv_threshold_4pl <- function(temp, survival_target,
                                       low, up, k,
                                       mid_int, mid_temp, temp_mean) {
  temp_c <- temp - temp_mean
  mid    <- mid_int + mid_temp * temp_c
  if (is.character(survival_target) && length(survival_target) == 1L &&
      survival_target == "relative") {
    return(10 ^ mid)
  }
  if (survival_target <= low || survival_target >= up)
    return(rep(NA_real_, length(temp)))
  log_term <- log((up - survival_target) / (survival_target - low))
  10 ^ (mid + log_term / k)
}

#' Survival corresponding to an accumulated dose
#'
#' Inverts the dose interpretation: at `dose = 1`, survival equals the
#' threshold by construction. The 4PL maps dose smoothly to survival on the
#' reference dose-response curve.
#'
#' @param dose Numeric vector of cumulative doses (in units where 1 dose =
#'             1 LT-dose at the chosen threshold).
#' @param low,up,k 4PL parameters at the reference (centring) temperature.
#' @param target_surv Either the literal probability defining "1 dose" or the
#'                    string `"relative"` for the `(low + up)/2` threshold.
#' @return Numeric vector of predicted survival probabilities.
#' @keywords internal
survival_from_dose <- function(dose, low, up, k, target_surv = "relative") {
  dose_use <- pmax(dose, 1e-12)
  if (is.character(target_surv) && length(target_surv) == 1L &&
      target_surv == "relative") {
    c_target <- 0
  } else {
    c_target <- log((up - target_surv) / (target_surv - low)) / k
  }
  low + (up - low) / (1 + exp(k * (log10(dose_use) + c_target)))
}

#' Predict heat injury and survival under a fluctuating temperature trace
#'
#' Propagates the model's full posterior through an Eulerian damage-
#' accumulation integration along the supplied temperature trace, returning
#' the posterior median and 95% credible band of:
#'
#' - **HI(t)** — cumulative heat injury, in percent of an LT_target_surv dose.
#'   When `HI(t) = 100`, the population has accumulated one full dose at the
#'   chosen survival threshold (default 50% mortality).
#' - **S(t)** — predicted survival fraction, mapped from the cumulative dose
#'   through the fitted 4PL.
#'
#' Optionally adds a temperature-dependent repair rate via
#' [repair_rate_schoolfield()].
#'
#' @param trace        Tibble with columns `time` (numeric time from start,
#'                     in `trace_unit`) and `temp` (°C), in time order.
#'                     Requires ≥ 2 rows.
#' @param workflow     Fitted `bayes_tls`.
#' @param target_surv  Threshold defining "1 dose". `"relative"` (default;
#'                     `(low + up)/2`), `"absolute"` (= 0.5), or a numeric in
#'                     `(0, 1)`. The default coincides with the classical LT50
#'                     when `low ≈ 0` and `up ≈ 1`; with sub-unit asymptotes
#'                     the relative threshold is the more biologically
#'                     meaningful anchor for dose accounting.
#' @param T_c          Optional damage-accumulation threshold (°C). When
#'                     supplied, the damage rate is forced to zero at
#'                     `temp <= T_c` (matches Equation 7 of the manuscript).
#'                     Default `NULL` lets the rate fall naturally with T.
#' @param trace_unit   Time unit of the trace's `time` column: one of
#'                     `"hours"` (default), `"minutes"`, `"seconds"`, `"days"`.
#'                     Reconciled internally with the model's fitted
#'                     `duration_unit`, so the result is correct for any
#'                     combination of model and trace time units.
#' @param ndraws       Posterior draws to use. Default 500.
#' @param repair       Logical. If `TRUE`, add Sharpe-Schoolfield repair.
#'                     Default `FALSE`.
#' @param repair_pars  Required when `repair = TRUE`. Named list with elements
#'                     `TA, TAL, TAH, TL, TH, TREF, r_ref` passed straight to
#'                     [repair_rate_schoolfield()]. `r_ref` should be in
#'                     "doses per hour" so it matches the damage-rate units.
#' @param repair_scales_with_survival Logical. If `TRUE` (default), repair
#'                     rate at each step is scaled by `survival / up` so dead
#'                     organisms do not contribute to repair.
#' @param irreversible_mortality Logical. If `TRUE` (default), survival can
#'                     only decrease over time — i.e. once the population's
#'                     predicted survival reaches a value, it cannot rebound
#'                     even if cumulative dose subsequently decreases.
#' @param save_draws   Logical. If `TRUE`, return the full per-draw
#'                     trajectories. Default `FALSE`.
#' @param seed         Optional integer seeding the posterior-draw subsample for
#'                     reproducibility. `NULL` (default) leaves the RNG untouched.
#' @param by           Optional moderator column(s) for per-group injury. `NULL`
#'                     (default) uses the fit's moderators; a single-condition fit
#'                     returns the ungrouped result. A grouped fit runs the dose
#'                     integral through each group's 4PL and `summary` gains the
#'                     moderator column(s).
#' @return A list with elements:
#'   - `summary`: tibble with `time`, `temp`, and posterior median + 95%
#'     CrI for `hi`, `survival`, and `mortality` at each time step (plus the
#'     moderator column(s) for a grouped fit).
#'   - `draws`: optional per-draw trajectories (when `save_draws = TRUE`).
#'   - `meta`: list of inputs used.
#' @examples
#' \dontrun{
#' scens <- make_temperature_scenarios()
#' hi    <- predict_heat_injury(scens$single_spike, wf)
#' hi$summary
#' }
#' @export
predict_heat_injury <- function(trace, workflow,
                                target_surv = "relative",
                                T_c         = NULL,
                                trace_unit  = "hours",
                                ndraws      = 500,
                                repair      = FALSE,
                                repair_pars = NULL,
                                repair_scales_with_survival = TRUE,
                                irreversible_mortality      = TRUE,
                                save_draws                  = FALSE,
                                seed                        = NULL,
                                by                          = NULL) {

  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  if (nrow(trace) < 2L)
    stop("trace must have at least 2 rows.", call. = FALSE)
  if (!all(c("time", "temp") %in% names(trace)))
    stop("trace must have `time` and `temp` columns.", call. = FALSE)
  if (anyNA(trace$time) || anyNA(trace$temp))
    stop("trace has NA in `time` or `temp`. Interpolate or drop missing rows ",
         "first -- an NA temperature would otherwise be silently counted as zero ",
         "heat injury, under-counting the dose.", call. = FALSE)
  if (repair && is.null(repair_pars))
    stop("Supply repair_pars when repair = TRUE.", call. = FALSE)

  # Normalise the threshold up front so we can pass either the literal
  # probability or the sentinel string "relative" down to the helpers.
  ts <- resolve_target_surv(target_surv)
  ts_arg <- if (ts$mode == "relative") "relative" else ts$prob

  trace     <- trace[order(trace$time), , drop = FALSE]
  n         <- nrow(trace)
  temp_mean <- workflow$meta$temp_mean

  # Reconcile BOTH time units to a common base (hours), so the integral is
  # correct no matter how the model was fit OR how the trace is labelled:
  #   - `tau` is returned in the model's `duration_unit`;
  #   - the trace step is in `trace_unit` (the `time` column).
  # We convert the damage rate to "doses per hour" and the step `dt` to hours.
  # Without this, e.g. a minutes-fitted model driven by an hours trace would
  # under-count the accumulated dose 60-fold (and the reverse over-counts).
  to_hours <- function(u) switch(as.character(u),
                   seconds = 1 / 3600, minutes = 1 / 60, hours = 1, days = 24,
                   stop(sprintf(
                     "predict_heat_injury(): unsupported time unit '%s'; expected one of seconds/minutes/hours/days.",
                     u), call. = FALSE))
  unit_h <- to_hours(workflow$meta$duration_unit %||% "hours")               # model -> hours
  # Per-interval widths (hours). Using diff() rather than reusing a single first
  # step is what makes irregular / gappy traces integrate correctly; dt_vec[j-1]
  # is the width of the interval ending at point j.
  dt_vec <- diff(trace$time) * to_hours(trace_unit)                          # trace -> hours

  by       <- tdt_resolve_by(workflow, by)
  pars_all <- extract_4pl_pars(workflow, by = by)

  # Per-draw Euler dose-accumulation integral for one set of 4PL parameter draws.
  integrate_pars <- function(pars) {
    pred_list <- vector("list", nrow(pars))
    for (i in seq_len(nrow(pars))) {
      tau <- time_to_surv_threshold_4pl(
        temp = trace$temp, survival_target = ts_arg,
        low = pars$low[i], up = pars$up[i], k = pars$k[i],
        mid_int = pars$mid_int[i], mid_temp = pars$mid_temp[i],
        temp_mean = temp_mean
      )
      dmg <- 1 / (tau * unit_h)               # doses per hour (unit-reconciled)
      dmg[!is.finite(dmg)] <- 0
      if (!is.null(T_c)) dmg[trace$temp <= T_c] <- 0

      rep_rate <- if (repair) {
        repair_rate_schoolfield(
          temp_celsius = trace$temp,
          TA = repair_pars$TA, TAL = repair_pars$TAL,
          TAH = repair_pars$TAH, TL = repair_pars$TL,
          TH = repair_pars$TH, TREF = repair_pars$TREF,
          r_ref = repair_pars$r_ref
        )
      } else {
        rep(0, n)
      }

      # Forward-Euler dose ODE: d(dose)/dt = damage(T) - repair(T) * scale. The
      # dose at the first time point is zero (no exposure has elapsed); each later
      # point adds the rate at the START of its interval times that interval's
      # width. (The previous loop credited a full step at t = 0 and reused the
      # first dt for every interval -- over-counting by one step and corrupting
      # irregular traces.) Kept consistent with planted_dose_from_trace().
      dose <- numeric(n); survival <- numeric(n)
      dose[1]     <- 0
      survival[1] <- survival_from_dose(0, low = pars$low[i], up = pars$up[i],
                                        k = pars$k[i], target_surv = ts_arg)
      for (j in seq_len(n)[-1]) {
        w     <- dt_vec[j - 1]
        rep_j <- rep_rate[j - 1] * w
        if (repair_scales_with_survival) rep_j <- rep_j * survival[j - 1] / pars$up[i]
        new_dose <- max(0, dose[j - 1] + dmg[j - 1] * w - rep_j)
        surv_raw <- survival_from_dose(
          new_dose, low = pars$low[i], up = pars$up[i], k = pars$k[i],
          target_surv = ts_arg
        )
        survival[j] <- if (irreversible_mortality) min(survival[j - 1], surv_raw) else surv_raw
        dose[j]     <- new_dose
      }
      pred_list[[i]] <- data.frame(
        .draw = pars$.draw[i], time = trace$time, temp = trace$temp,
        dose = dose, hi = dose * 100, survival = survival, mortality = 1 - survival
      )
    }
    dplyr::bind_rows(pred_list)
  }

  if (!is.null(seed)) set.seed(seed)   # reproducible posterior-draw subsample
  if (is.null(by)) {
    draws <- integrate_pars(dplyr::slice_sample(pars_all, n = min(ndraws, nrow(pars_all))))
    n_used <- min(ndraws, nrow(pars_all))
  } else {
    groups <- unique(pars_all[, by, drop = FALSE])
    draws <- dplyr::bind_rows(lapply(seq_len(nrow(groups)), function(gi) {
      gp <- dplyr::inner_join(pars_all, groups[gi, , drop = FALSE], by = by)
      cbind(groups[gi, , drop = FALSE],
            integrate_pars(dplyr::slice_sample(gp, n = min(ndraws, nrow(gp)))),
            row.names = NULL)
    }))
    n_used <- ndraws
  }

  q_lower <- function(x) stats::quantile(x, 0.025, na.rm = TRUE)
  q_upper <- function(x) stats::quantile(x, 0.975, na.rm = TRUE)
  summary <- draws |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(by, "time", "temp")))) |>
    dplyr::summarise(
      hi_median   = stats::median(hi,        na.rm = TRUE),
      hi_lower    = q_lower(hi),
      hi_upper    = q_upper(hi),
      surv_median = stats::median(survival,  na.rm = TRUE),
      surv_lower  = q_lower(survival),
      surv_upper  = q_upper(survival),
      mort_median = stats::median(mortality, na.rm = TRUE),
      mort_lower  = q_lower(mortality),
      mort_upper  = q_upper(mortality),
      .groups     = "drop"
    )

  out <- list(
    summary = summary,
    meta    = list(target_surv = ts$label, target_mode = ts$mode,
                   target_prob = ts$prob, T_c = T_c, repair = repair,
                   repair_pars = repair_pars, ndraws = n_used, by = by)
  )
  if (save_draws) out$draws <- draws
  out
}
