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

#' Per-draw natural-scale 4PL parameters from a fitted workflow
#'
#' Faithfully extracts ALL FOUR natural-scale 4PL sub-parameters
#' (`low`, `up`, `k`, `mid`) per posterior draw, at one or more assay
#' temperatures and per moderator group, retaining whatever temperature and/or
#' group structure each sub-parameter was fit with (`~ 1`, `~ temp_c`,
#' `~ temp_c * group`, random effects, …). It surfaces the shared TDT engine
#' (`tls_eval_subpars()` over a `tls_build_grid()` temperature × moderator
#' grid via `brms::posterior_linpred(nlpar=)`), so the result is
#' parameterisation- and coding-agnostic — no coefficient-name parsing.
#'
#' `mid` is the natural-scale midpoint sub-parameter (the log10 time at the
#' per-draw relative threshold) evaluated at each temperature; on an absolute
#' fit it is still the bare `mid` nlpar (the asymmetry correction is applied
#' downstream by the threshold-specific readers, not folded in here).
#'
#' This replaces the earlier constant-shape extractor that discarded the
#' temperature slopes on `low`/`up`/`k`. The constant-in-T view the
#' heat-injury integral needs (`low, up, k` at the centring temperature plus
#' a linear `mid_int`/`mid_temp`) now lives in the internal helper
#' [hi_pars()], which [predict_heat_injury()] calls for `shape = "constant"`.
#'
#' @param workflow Fitted `bayes_tls`.
#' @param temps Numeric vector of assay temperatures (°C, uncentred) to
#'   evaluate the curve parameters at. `NULL` (default) uses the fit's observed
#'   unique assay temperatures.
#' @param by Optional moderator column(s) for per-group parameters. `NULL`
#'   (default) uses the fit's moderators (all levels); a single-condition fit
#'   returns the ungrouped tibble. A grouped fit prepends the moderator
#'   column(s).
#' @param re_formula Passed to [brms::posterior_linpred()]. `NA` (default)
#'   marginalises out group-level (random) effects so the result is a
#'   population-level curve; `NULL` conditions on the fitted random effects.
#' @param ndraws Posterior draws to use; `NULL` (default) uses the full
#'   posterior. Pass an integer to subsample (the caller is responsible for
#'   `set.seed()` if a reproducible subsample is wanted).
#' @return A long tibble with `(.draw, temp, low, up, k, mid)` columns (plus
#'         the moderator column(s) for a grouped fit, immediately after
#'         `.draw`), one row per draw × temperature × group, filtered to rows
#'         producing valid parameter values (`k > 0`, `up > low`, all finite).
#'         `temp` is the uncentred assay temperature (°C).
#' @seealso [hi_pars()] for the constant-shape view used by the heat-injury
#'   integral; [get_4pl_est()] which wraps this for the `"draws"` view.
#' @examples
#' \dontrun{
#' pars <- extract_4pl_pars(wf)                 # observed temps, fit's groups
#' extract_4pl_pars(wf, temps = c(30, 33, 36))  # at chosen temperatures
#' head(pars)
#' }
#' @export
extract_4pl_pars <- function(workflow, temps = NULL, by = NULL,
                             re_formula = NA, ndraws = NULL) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  fit       <- get_brmsfit(workflow)
  by        <- tdt_resolve_by(workflow, by)
  temp_mean <- workflow$meta$temp_mean %||% 0

  # Default to the fit's observed unique assay temperatures (uncentred °C).
  if (is.null(temps)) {
    temps <- sort(unique(fit$data$temp_c)) + temp_mean
  }
  temp_grid <- temps - temp_mean                       # centre for the engine

  draw_ids <- if (is.null(ndraws)) NULL else tls_draw_ids(fit, ndraws)
  nd <- tls_build_grid(fit$data, by = by, temp = "temp_c", temp_grid = temp_grid)
  sp <- tls_eval_subpars(fit, nd, workflow$meta$bounds, mode = "relative",
                         re_formula = re_formula, draw_ids = draw_ids)

  ndr <- nrow(sp$mid)
  per_group <- lapply(unique(nd$.grp), function(g) {
    gi <- which(nd$.grp == g)
    # One block per (draw, temperature) cell in this group, in temp order.
    blocks <- lapply(gi, function(col) {
      tibble::tibble(
        .draw = seq_len(ndr),
        temp  = nd$temp_c[col] + temp_mean,
        low   = sp$low[, col], up = sp$up[, col], k = sp$k[, col],
        mid   = sp$mid[, col]
      )
    })
    out <- dplyr::bind_rows(blocks)
    out <- dplyr::filter(out,
                         is.finite(low) & is.finite(up) & is.finite(k) &
                         is.finite(mid) & k > 0 & up > low)
    if (!is.null(by)) {
      gcols <- nd[gi[1], by, drop = FALSE]
      out <- cbind(gcols[rep(1, nrow(out)), , drop = FALSE], out, row.names = NULL)
    }
    out
  })
  tibble::as_tibble(dplyr::bind_rows(per_group))
}

#' Constant-shape 4PL parameters for the heat-injury integral (internal)
#'
#' The classical heat-injury integral is evaluated under the assumption that
#' the asymptotes (`low`, `up`) and steepness (`k`) are constant in
#' temperature and only the midpoint shifts. This helper returns exactly that
#' view, per posterior draw and per moderator group: `low`, `up`, `k` read at
#' the centring temperature (`temp_c = 0`, i.e. `T_bar`), plus a linear
#' midpoint decomposed into intercept (`mid_int`, at `T_bar`) and slope
#' (`mid_temp`, per °C). This is the per-draw input
#' [predict_heat_injury()] uses for `shape = "constant"`.
#'
#' Reads `posterior_linpred(nlpar=)` at `temp_c = 0` (asymptotes + midpoint
#' intercept) and `temp_c = 1` (for the linear midpoint slope), so it is
#' parameterisation- and coding-agnostic. If the temperature slopes on
#' `low`/`up`/`k` are shrunk near zero by the data this introduces no bias;
#' otherwise it is the same approximation the classical HI framework makes
#' (see `shape = "varying"` for the slope-aware alternative).
#'
#' @param workflow Fitted `bayes_tls`.
#' @param by Resolved moderator column(s) (already through `tdt_resolve_by()`),
#'   or `NULL` for a single-condition fit.
#' @return A tibble with `(.draw, low, up, k, mid_int, mid_temp)` (plus the
#'   moderator column(s) for a grouped fit), filtered to valid draws.
#' @keywords internal
hi_pars <- function(workflow, by = NULL) {
  fit <- get_brmsfit(workflow)
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

#' Temperature-local 4PL curves per posterior draw (internal, for shape="varying")
#'
#' Evaluates the full natural-scale 4PL `low(T)`, `up(T)`, `k(T)`, `mid(T)` per
#' posterior draw at a set of unique temperatures, retaining every sub-parameter's
#' temperature/group structure (the slopes the constant-shape [hi_pars()] view
#' drops). Returns per-draw × per-temperature matrices, keyed for fast lookup by
#' the shape-varying heat-injury integrator: each trace point's temperature maps
#' to a column, each posterior draw to a row.
#'
#' @param workflow Fitted `bayes_tls`.
#' @param by Resolved moderator column(s), or `NULL` for a single-condition fit.
#' @param temps_unique Numeric vector of unique (uncentred °C) temperatures to
#'   evaluate the curves at -- typically `unique(trace$temp)`.
#' @return A named list, one element per moderator group (named by the group
#'   label, `"all"` when ungrouped). Each element is a list with:
#'   `temps` (the input temperatures, in evaluation order), `.draw` (the
#'   posterior draw indices kept after filtering invalid draws), and the
#'   matrices `low`, `up`, `k`, `mid`, each `[length(.draw) × length(temps)]`.
#'   A draw is kept only if ALL its cells are valid (finite, `k > 0`,
#'   `up > low`) at every temperature, so the integrator never hits a partial
#'   draw.
#' @keywords internal
hi_local_curves <- function(workflow, by = NULL, temps_unique) {
  fit       <- get_brmsfit(workflow)
  temp_mean <- workflow$meta$temp_mean %||% 0
  temp_grid <- temps_unique - temp_mean
  nd <- tls_build_grid(fit$data, by = by, temp = "temp_c", temp_grid = temp_grid)
  sp <- tls_eval_subpars(fit, nd, workflow$meta$bounds, mode = "relative")
  ndr <- nrow(sp$mid)

  out <- lapply(unique(nd$.grp), function(g) {
    gi <- which(nd$.grp == g)
    # Columns of nd for this group, ordered to match temps_unique.
    ord <- gi[order(match(round(nd$temp_c[gi] + temp_mean, 6),
                          round(temps_unique, 6)))]
    low <- sp$low[, ord, drop = FALSE]; up <- sp$up[, ord, drop = FALSE]
    k   <- sp$k[, ord, drop = FALSE];  mid <- sp$mid[, ord, drop = FALSE]
    # Keep only draws valid at EVERY temperature (so no partial-draw traces).
    ok <- apply(is.finite(low) & is.finite(up) & is.finite(k) & is.finite(mid) &
                  (k > 0) & (up > low), 1, all)
    list(temps = temps_unique, .draw = which(ok), grp = g,
         low = low[ok, , drop = FALSE], up = up[ok, , drop = FALSE],
         k = k[ok, , drop = FALSE],   mid = mid[ok, , drop = FALSE],
         by_row = if (!is.null(by)) nd[gi[1], by, drop = FALSE] else NULL)
  })
  names(out) <- vapply(out, function(x) x$grp, character(1))
  out
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
#' @return Numeric vector of predicted survival probabilities. `NA` for every
#'         element when a numeric `target_surv` lies outside `(low, up)` and is
#'         therefore unreachable on this draw's 4PL (mirrors
#'         [time_to_surv_threshold_4pl()], so the draw drops cleanly rather than
#'         producing `log()`-of-negative `NaN`s that would poison the trace).
#' @keywords internal
survival_from_dose <- function(dose, low, up, k, target_surv = "relative") {
  dose_use <- pmax(dose, 1e-12)
  if (is.character(target_surv) && length(target_surv) == 1L &&
      target_surv == "relative") {
    c_target <- 0
  } else {
    # Unreachable target on this draw -- guard like time_to_surv_threshold_4pl()
    # so the whole draw drops as NA instead of evaluating log() of a negative
    # number (which would yield NaN and silently bias the medians upward).
    if (target_surv <= low || target_surv >= up)
      return(rep(NA_real_, length(dose)))
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
#'                     Requires >= 2 rows.
#' @param workflow     Fitted `bayes_tls`.
#' @param target_surv  Threshold defining "1 dose". `"relative"` (default;
#'                     `(low + up)/2`), `"absolute"` (= 0.5), or a numeric in
#'                     `(0, 1)`. The default coincides with the classical LT50
#'                     when `low` is about 0 and `up` about 1; with sub-unit asymptotes
#'                     the relative threshold is the more biologically
#'                     meaningful anchor for dose accounting.
#' @param T_c          Optional damage-accumulation threshold (°C). When
#'                     supplied, the damage rate is forced to zero at
#'                     `temp <= T_c` (matches the heat-injury integral,
#'                     Equation 8 of the manuscript).
#'                     Default `NULL` lets the rate fall naturally with T.
#' @param trace_unit   Time unit of the trace's `time` column: one of
#'                     `"hours"` (default), `"minutes"`, `"seconds"`, `"days"`.
#'                     Reconciled internally with the model's fitted
#'                     `duration_unit`, so the result is correct for any
#'                     combination of model and trace time units.
#' @param ndraws       Posterior draws to use; `NULL` (default) uses the full
#'                     posterior, keeping the integral consistent with the model.
#'                     Pass an integer to subsample -- useful for speed, since the
#'                     full posterior over a long temperature trace can be slow.
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
#' @param shape        How the 4PL asymptotes/steepness enter the dose integral.
#'                     `"constant"` (default) holds `low`, `up`, `k` at the
#'                     centring temperature `T_bar` and lets only the midpoint
#'                     shift with temperature -- the classical constant-shape
#'                     heat-injury assumption (via [hi_pars()]). `"varying"`
#'                     feeds the temperature-local `low(T)`, `up(T)`, `k(T)`
#'                     (and `mid(T)`) into both the damage rate and the
#'                     dose→survival mapping (via [extract_4pl_pars()]). In
#'                     `"varying"` mode survival is accumulated as a monotone
#'                     state by an incremental decrement: each interval's
#'                     mortality increment is the *local* curve's drop in
#'                     survival over the dose added that interval
#'                     (`g(D_j; shape_j) - g(D_{j-1}; shape_j)`), so the local
#'                     shape sets the RATE of survival loss without re-mapping
#'                     the global cumulative dose through a step-local curve
#'                     (which would make survival jump when temperature changes
#'                     with no added exposure). When the shape is flat in
#'                     temperature the two modes coincide to numerical tolerance.
#'                     See `vignette`-free note `notes/2026-06-26-shape-varying-heat-injury.qmd`.
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
                                ndraws      = NULL,
                                repair      = FALSE,
                                repair_pars = NULL,
                                repair_scales_with_survival = TRUE,
                                irreversible_mortality      = TRUE,
                                save_draws                  = FALSE,
                                seed                        = NULL,
                                by                          = NULL,
                                shape = c("constant", "varying")) {

  shape <- match.arg(shape)
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
  pars_all <- hi_pars(workflow, by = by)
  # ndraws = NULL (default) uses the FULL posterior -- keeps the integral
  # consistent with the model and the other derivation functions; pass an
  # integer to subsample for speed on long traces.
  if (is.null(ndraws)) ndraws <- nrow(pars_all)

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

  # ----- shape = "varying": temperature-LOCAL low/up/k/mid in the integral -----
  # The local curves are evaluated once per group at the trace's unique temps;
  # `curve$.draw` indexes the full posterior, so a sampled draw is matched by its
  # .draw id. `tcol` maps each trace point to its temperature column.
  temps_unique <- sort(unique(trace$temp))
  tcol         <- match(trace$temp, temps_unique)
  curves_by_grp <- if (shape == "varying")
    hi_local_curves(workflow, by = by, temps_unique = temps_unique) else NULL

  # Survival at cumulative dose D on a 4PL curve (low, up, k) anchored so that
  # D = 1 hits the chosen threshold (relative midpoint, or absolute p). Same map
  # as survival_from_dose(), vectorised over per-step (D, low, up, k). Used to
  # form the LOCAL per-interval survival decrement g(D_j) - g(D_{j-1}) (both
  # endpoints on the SAME local curve, so a temperature change adds no jump).
  surv_at_dose <- function(D, low, up, k) {
    D <- pmax(D, 1e-12)
    if (is.character(ts_arg) && length(ts_arg) == 1L && ts_arg == "relative") {
      c_t <- 0
    } else {
      ratio <- (up - ts_arg) / (ts_arg - low)
      c_t   <- ifelse(is.finite(ratio) & ratio > 0, log(ratio) / k, NA_real_)
    }
    low + (up - low) / (1 + exp(k * (log10(D) + c_t)))
  }

  integrate_pars_varying <- function(pars, curve) {
    # Row in the local-curve matrices for each sampled draw (matched by .draw id).
    rmap <- match(pars$.draw, curve$.draw)
    pred_list <- vector("list", nrow(pars))
    for (i in seq_len(nrow(pars))) {
      r <- rmap[i]
      if (is.na(r)) next                      # draw dropped as invalid somewhere
      lowT <- curve$low[r, tcol]; upT <- curve$up[r, tcol]
      kT   <- curve$k[r, tcol];   midT <- curve$mid[r, tcol]

      # Damage rate from the LOCAL curve: tau(T) via the natural-scale mid(T)
      # (mid_int = mid(T), mid_temp = 0 -> mid evaluated exactly at each point),
      # with the local low/up/k feeding the absolute-threshold inversion.
      tau <- vapply(seq_len(n), function(j) time_to_surv_threshold_4pl(
        temp = trace$temp[j], survival_target = ts_arg,
        low = lowT[j], up = upT[j], k = kT[j],
        mid_int = midT[j], mid_temp = 0, temp_mean = temp_mean), numeric(1))
      dmg <- 1 / (tau * unit_h)
      dmg[!is.finite(dmg)] <- 0
      if (!is.null(T_c)) dmg[trace$temp <= T_c] <- 0

      rep_rate <- if (repair) {
        repair_rate_schoolfield(
          temp_celsius = trace$temp,
          TA = repair_pars$TA, TAL = repair_pars$TAL,
          TAH = repair_pars$TAH, TL = repair_pars$TL,
          TH = repair_pars$TH, TREF = repair_pars$TREF,
          r_ref = repair_pars$r_ref)
      } else rep(0, n)

      dose <- numeric(n); survival <- numeric(n)
      dose[1]     <- 0
      survival[1] <- surv_at_dose(1e-12, lowT[1], upT[1], kT[1])   # = up at D->0
      for (j in seq_len(n)[-1]) {
        w     <- dt_vec[j - 1]
        rep_j <- rep_rate[j - 1] * w
        if (repair_scales_with_survival) rep_j <- rep_j * survival[j - 1] / upT[j]
        new_dose <- max(0, dose[j - 1] + dmg[j - 1] * w - rep_j)
        # Incremental survival decrement from the LOCAL curve over [D_{j-1}, D_j]:
        # both doses evaluated on the SAME step-j curve, so the shape change
        # contributes no discontinuity; the local slope-in-dose sets the rate of
        # survival loss. When low/up/k are flat in T this telescopes EXACTLY to
        # the constant-shape closed form g(D_n).
        s_prev_local <- surv_at_dose(dose[j - 1], lowT[j], upT[j], kT[j])
        s_new_local  <- surv_at_dose(new_dose,    lowT[j], upT[j], kT[j])
        dec <- s_new_local - s_prev_local                  # <= 0 (dose increases)
        if (!is.finite(dec)) dec <- 0
        surv_raw    <- survival[j - 1] + dec
        survival[j] <- if (irreversible_mortality) min(survival[j - 1], surv_raw) else surv_raw
        survival[j] <- max(survival[j], 0)                 # respect the floor
        dose[j]     <- new_dose
      }
      pred_list[[i]] <- data.frame(
        .draw = pars$.draw[i], time = trace$time, temp = trace$temp,
        dose = dose, hi = dose * 100, survival = survival, mortality = 1 - survival
      )
    }
    dplyr::bind_rows(pred_list)
  }

  # Dispatch: pick the per-draw integrator for the requested shape. The grouped
  # path selects the matching group's local curves by label.
  run_integral <- function(pars, grp_label) {
    if (shape == "constant") integrate_pars(pars)
    else integrate_pars_varying(pars, curves_by_grp[[grp_label]])
  }

  if (!is.null(seed)) set.seed(seed)   # reproducible posterior-draw subsample
  if (is.null(by)) {
    draws <- run_integral(
      dplyr::slice_sample(pars_all, n = min(ndraws, nrow(pars_all))), "all")
    n_used <- min(ndraws, nrow(pars_all))
  } else {
    groups <- unique(pars_all[, by, drop = FALSE])
    per_group_n <- integer(nrow(groups))
    draws <- dplyr::bind_rows(lapply(seq_len(nrow(groups)), function(gi) {
      gp     <- dplyr::inner_join(pars_all, groups[gi, , drop = FALSE], by = by)
      n_gi   <- min(ndraws, nrow(gp))
      per_group_n[gi] <<- n_gi
      # Group label, built EXACTLY as tls_build_grid() builds its `.grp`
      # (do.call(paste, c(<by columns>, sep = " / "))) so the varying-mode local
      # curves -- keyed by that same `.grp` -- are matched for this group.
      grp_label <- do.call(paste, c(groups[gi, by, drop = FALSE], sep = " / "))
      cbind(groups[gi, , drop = FALSE],
            run_integral(dplyr::slice_sample(gp, n = n_gi), grp_label),
            row.names = NULL)
    }))
    # Report the actual clamped per-group count integrated, not the raw request.
    n_used <- if (length(unique(per_group_n)) == 1L) per_group_n[1]
              else sum(per_group_n)
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
                   repair_pars = repair_pars, ndraws = n_used, by = by,
                   shape = shape)
  )
  if (save_draws) out$draws <- draws
  out
}
