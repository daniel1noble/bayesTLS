# Heat injury and predicted survival under a fluctuating temperature trace.
# Carries the model's full posterior through to HI(t) and S(t) trajectories
# with credible intervals. The algorithm:
#
#   1. Extract per-draw 4PL parameters (low, up, k, mid_int, mid_temp).
#   2. At each time step, compute LT_target(T) analytically from the 4PL,
#      giving damage rate = 1 / LT_target.
#   3. Optionally add a temperature-dependent repair rate (Sharpe-Schoolfield).
#   4. Trapezoidally integrate (damage - repair) across the trace.
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
#' @return A tibble with `(.draw, low, up, k, mid_int, mid_temp)` columns,
#'         filtered to draws producing valid parameter values.
#' @examples
#' \dontrun{
#' pars <- extract_4pl_pars(wf)
#' head(pars)
#' }
#' @export
extract_4pl_pars <- function(workflow) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  d  <- posterior::as_draws_df(workflow$fit) |> as.data.frame()
  b  <- workflow$meta$bounds

  # Direct fits carry CTmaxdev/logz, not b_mid_*; a grouped direct fit has no
  # *_Intercept coefficients at all. Guard before touching any Intercept so the
  # message is clear rather than a downstream plogis(NULL) error.
  direct <- identical(workflow$meta$parameterization, "direct")
  if (direct && !all(c("b_CTmaxdev_Intercept", "b_logz_Intercept") %in% names(d)))
    stop("extract_4pl_pars(): single-group direct fits only; for grouped ",
         "direct fits use tls() / posterior_linpred-based helpers.",
         call. = FALSE)

  low <- b$low_min + stats::plogis(d$b_lowraw_Intercept) * b$low_w
  up  <- b$up_min  + stats::plogis(d$b_upraw_Intercept)  * b$up_w
  k   <- exp(d$b_logk_Intercept)

  # Linear midpoint coefficients (mid = mid_int + mid_temp * temp_c). Midpoint
  # parameterisation reads them directly. Direct parameterisation has no `mid`
  # coefficient, so reconstruct the *relative* midpoint backbone from
  # CTmaxdev/logz: mid(T) = log10_tref - (temp_c - CTmaxdev)/exp(logz), giving
  # slope -1/exp(logz) and intercept log10_tref + CTmaxdev/exp(logz). Under the
  # heat-injury constant-shape assumption (low/up/k read at the Intercept), an
  # absolute fit's relative midpoint is that backbone shifted by a constant
  # C = log((up - 0.5)/(0.5 - low))/k, so it stays linear here.
  if (direct) {
    z        <- exp(d$b_logz_Intercept)
    l10      <- workflow$meta$log10_tref %||% 0
    mid_temp <- -1 / z
    mid_int  <- l10 + d$b_CTmaxdev_Intercept / z
    if (identical(workflow$meta$threshold %||% "relative", "absolute"))
      mid_int <- mid_int - log((up - 0.5) / (0.5 - low)) / k
  } else {
    mid_int  <- d$b_mid_Intercept
    mid_temp <- d$b_mid_temp_c
  }

  out <- tibble::tibble(
    .draw    = seq_len(nrow(d)),
    low      = low,
    up       = up,
    k        = k,
    mid_int  = mid_int,
    mid_temp = mid_temp
  )

  dplyr::filter(out,
                is.finite(low) & is.finite(up) & is.finite(k) &
                is.finite(mid_int) & is.finite(mid_temp) &
                k > 0 & up > low)
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
#' @return A list with elements:
#'   - `summary`: tibble with `time`, `temp`, and posterior median + 95%
#'     CrI for `hi`, `survival`, and `mortality` at each time step.
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
                                save_draws                  = FALSE) {

  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  if (nrow(trace) < 2L)
    stop("trace must have at least 2 rows.", call. = FALSE)
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
  dt     <- (if (n >= 2L) diff(trace$time)[1] else 1) * to_hours(trace_unit)  # trace -> hours

  pars <- extract_4pl_pars(workflow)
  pars <- dplyr::slice_sample(pars, n = min(ndraws, nrow(pars)))

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

    # Euler integration matching planted_dose_from_trace() and supplement
    # Equation 7: rate at point j applies for dt forward; cumulative dose at
    # point i is the sum of (rate × dt) for j = 1..i.
    dose      <- numeric(n)
    survival  <- numeric(n)
    prev_dose <- 0
    prev_surv <- pars$up[i]

    for (j in seq_len(n)) {
      rep_j <- rep_rate[j] * dt
      if (repair_scales_with_survival) rep_j <- rep_j * prev_surv / pars$up[i]

      new_dose <- max(0, prev_dose + dmg[j] * dt - rep_j)
      surv_raw <- survival_from_dose(
        new_dose, low = pars$low[i], up = pars$up[i], k = pars$k[i],
        target_surv = ts_arg
      )
      surv_new <- if (irreversible_mortality) min(prev_surv, surv_raw) else surv_raw

      dose[j]     <- new_dose
      survival[j] <- surv_new
      prev_dose   <- new_dose
      prev_surv   <- surv_new
    }

    pred_list[[i]] <- data.frame(
      .draw     = pars$.draw[i],
      time    = trace$time,
      temp      = trace$temp,
      dose      = dose,
      hi        = dose * 100,        # % of LT_{target_surv} dose
      survival  = survival,
      mortality = 1 - survival
    )
  }

  draws <- dplyr::bind_rows(pred_list)

  q_lower <- function(x) stats::quantile(x, 0.025, na.rm = TRUE)
  q_upper <- function(x) stats::quantile(x, 0.975, na.rm = TRUE)
  summary <- draws |>
    dplyr::group_by(time, temp) |>
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
                   repair_pars = repair_pars, ndraws = nrow(pars))
  )
  if (save_draws) out$draws <- draws
  out
}
