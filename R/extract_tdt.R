# Classical TDT quantities (z, CTmax, T_crit) derived from a fitted 4PL via
# numerical inversion of the posterior survival surface. Three primitives plus
# a bundling wrapper.

#' Posterior LT_x curve: time to reach a given survival target at each temperature
#'
#' For each posterior draw and each temperature in `temp_grid`, finds the
#' duration at which population-level survival crosses `target_surv`. The
#' inversion is numerical — predict survival on a dense duration grid, then
#' `approx()` through `target_surv`.
#'
#' This is the **horizontal** read of the survival surface: fix a survival
#' threshold, read off the time required to reach it at each temperature.
#'
#' @param workflow         Fitted `tdt_4pl_workflow`.
#' @param temp_grid        Numeric vector of temperatures (°C).
#' @param duration_grid    Numeric vector of durations along which to search.
#'                         Default: 350 log-spaced values spanning 0.2× to 5×
#'                         the training data's duration range.
#' @param target_surv      Survival probability to invert at. Default 0.5.
#' @param ndraws           Posterior draws to use. Default 1000.
#' @param probs            Quantile probabilities for the summary. Default
#'                         `c(0.025, 0.5, 0.975)`.
#' @param time_multiplier  Multiplier from model time units to output time
#'                         units (e.g. 60 for hours → min). Default 60.
#' @param output_time_unit Label for the output time unit. Default `"min"`.
#' @return A list with `draws` (per-draw threshold durations), `summary`
#'         (quantile summary by temperature), `target_surv`, `time_multiplier`,
#'         `output_time_unit`.
#' @export
derive_ltx_curve <- function(workflow,
                             temp_grid,
                             duration_grid    = NULL,
                             target_surv      = 0.5,
                             ndraws           = 1000,
                             probs            = c(0.025, 0.5, 0.975),
                             time_multiplier  = 60,
                             output_time_unit = "min") {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  if (is.null(duration_grid)) {
    drange <- range(workflow$data$duration, na.rm = TRUE)
    duration_grid <- 10 ^ seq(log10(drange[1] / 5),
                              log10(drange[2] * 5),
                              length.out = 350)
  }

  nd   <- new_tdt_grid(workflow, temps = temp_grid, durations = duration_grid)
  pred <- posterior_linpred_tdt(workflow, nd, ndraws = ndraws, re_formula = NA)

  draw_list <- vector("list", length(temp_grid))
  for (i in seq_along(temp_grid)) {
    t_i <- temp_grid[i]
    idx <- nd$temp == t_i
    thr <- threshold_x_by_draw(pred_mat = pred[, idx, drop = FALSE],
                               x        = nd$duration[idx],
                               target   = target_surv)
    draw_list[[i]] <- data.frame(
      .draw            = seq_along(thr),
      temp             = t_i,
      target_surv      = target_surv,
      duration_model   = thr,
      duration_out     = thr * time_multiplier
    )
  }

  draws <- dplyr::bind_rows(draw_list) |>
    dplyr::filter(is.finite(duration_model), duration_model > 0)

  summary <- draws |>
    dplyr::group_by(target_surv, temp) |>
    dplyr::summarise(
      duration_lower  = stats::quantile(duration_out, probs[1], na.rm = TRUE),
      duration_median = stats::quantile(duration_out, probs[2], na.rm = TRUE),
      duration_upper  = stats::quantile(duration_out, probs[3], na.rm = TRUE),
      .groups = "drop"
    )

  list(draws            = draws,
       summary          = summary,
       target_surv      = target_surv,
       time_multiplier  = time_multiplier,
       output_time_unit = output_time_unit)
}

#' Temperature at which survival equals a target after a fixed exposure
#'
#' The **vertical** read of the survival surface: fix an exposure duration,
#' invert numerically over temperature to find where survival crosses
#' `target_surv`. Returns one temperature per posterior draw.
#'
#' This is the primitive used by [extract_tdt()] to derive CTmax (at
#' `target_surv = 0.5`) and T_crit (at `target_surv = 1 - TC_thresh`).
#'
#' @param workflow         Fitted `tdt_4pl_workflow`.
#' @param exposure_duration Numeric scalar — the fixed duration (model units).
#' @param temp_grid        Numeric vector of temperatures to search over.
#' @param target_surv      Survival probability to invert at. Default 0.5.
#' @param ndraws           Posterior draws. Default 1000.
#' @param probs            Quantile probabilities. Default `c(0.025, 0.5, 0.975)`.
#' @return A list with `draws` (per-draw threshold temperatures), `summary`
#'         (quantile summary), `exposure_duration`, `target_surv`.
#' @export
derive_temperature_for_duration <- function(workflow,
                                            exposure_duration,
                                            temp_grid,
                                            target_surv = 0.5,
                                            ndraws      = 1000,
                                            probs       = c(0.025, 0.5, 0.975)) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  nd   <- new_tdt_grid(workflow, temps = temp_grid,
                       durations = exposure_duration)
  pred <- posterior_linpred_tdt(workflow, nd, ndraws = ndraws, re_formula = NA)

  thr <- threshold_x_by_draw(pred_mat = pred,
                             x        = nd$temp,
                             target   = target_surv)

  draws <- tibble::tibble(
    .draw       = seq_along(thr),
    target_surv = target_surv,
    temp        = thr
  ) |>
    dplyr::filter(is.finite(temp))

  summary <- draws |>
    dplyr::summarise(
      target_surv = target_surv[1],
      temp_lower  = stats::quantile(temp, probs[1], na.rm = TRUE),
      temp_median = stats::quantile(temp, probs[2], na.rm = TRUE),
      temp_upper  = stats::quantile(temp, probs[3], na.rm = TRUE)
    )

  list(draws             = draws,
       summary           = summary,
       exposure_duration = exposure_duration,
       target_surv       = target_surv)
}

#' Per-draw z and CTmax from an LT_x curve via per-draw log-linear regression
#'
#' Bayesian analogue of fitting the classical TDT line. For each posterior
#' draw, regress `log10(LT_x duration)` on temperature and read off the
#' classical TDT quantities:
#'
#' - `z = -1 / slope` — thermal sensitivity, in °C per decade of time
#' - `CTmax = (log10(t_ref) - intercept) / slope` — temperature at which
#'   LT_x equals the reference time.
#'
#' @param ltx_curve  Output of [derive_ltx_curve()].
#' @param t_ref      Reference time (in `output_time_unit` of `ltx_curve`) for
#'                   computing CTmax. Default 60 (i.e. CTmax_1hr in min units).
#' @param min_points Minimum number of temperatures a draw must have to be
#'                   included. Default 5.
#' @return A list with `draws` (per-draw slope, z, CTmax, R²), `summary`,
#'         `t_ref`.
#' @export
derive_tdt_parameters <- function(ltx_curve,
                                  t_ref      = 60,
                                  min_points = 5) {
  draws <- ltx_curve$draws |>
    dplyr::filter(is.finite(duration_out), duration_out > 0) |>
    dplyr::mutate(log10_duration = log10(duration_out))

  params <- draws |>
    dplyr::group_by(target_surv, .draw) |>
    dplyr::filter(dplyr::n() >= min_points) |>
    dplyr::group_modify(function(d, key) {
      fit <- stats::lm(log10_duration ~ temp, data = d)
      co  <- stats::coef(fit)
      data.frame(
        intercept = unname(co[1]),
        slope_T   = unname(co[2]),
        z         = -1 / unname(co[2]),
        CTmax     = (log10(t_ref) - unname(co[1])) / unname(co[2]),
        r_squared = summary(fit)$r.squared
      )
    }) |>
    dplyr::ungroup()

  summary <- params |>
    dplyr::group_by(target_surv) |>
    dplyr::summarise(
      slope_median = stats::median(slope_T, na.rm = TRUE),
      slope_lower  = stats::quantile(slope_T, 0.025, na.rm = TRUE),
      slope_upper  = stats::quantile(slope_T, 0.975, na.rm = TRUE),
      z_median     = stats::median(z, na.rm = TRUE),
      z_lower      = stats::quantile(z, 0.025, na.rm = TRUE),
      z_upper      = stats::quantile(z, 0.975, na.rm = TRUE),
      CTmax_median = stats::median(CTmax, na.rm = TRUE),
      CTmax_lower  = stats::quantile(CTmax, 0.025, na.rm = TRUE),
      CTmax_upper  = stats::quantile(CTmax, 0.975, na.rm = TRUE),
      r2_median    = stats::median(r_squared, na.rm = TRUE),
      .groups      = "drop"
    )

  list(draws = params, summary = summary, t_ref = t_ref)
}

#' Extract classical TDT quantities from a fitted 4PL: z, CTmax, T_crit
#'
#' Bundles the three classical TDT summaries with full posterior uncertainty:
#'
#' - **z** — thermal sensitivity, derived by per-draw log-linear regression of
#'   `log10(LT50)` on temperature (see [derive_tdt_parameters()]).
#' - **CTmax** — temperature at which survival = 0.5 after `t_ref` exposure;
#'   the temperature where the LT_50 curve crosses `t_ref`. Computed by
#'   inverting the 4PL surface at fixed exposure duration.
#' - **T_crit** — the temperature below which the TDT-line damage rate falls
#'   beneath an empirically motivated floor `r*` (% HI per hour). For each
#'   posterior draw, `T_crit = CTmax + z * log10(r* / 100)`, with `r*` drawn
#'   uniformly on the log10 scale across `TC_rate_range`. The pooled posterior
#'   thus carries both parameter uncertainty (in CTmax and z) and operational
#'   uncertainty in the choice of damage-rate floor. The default range
#'   `c(0.1, 1)` % HI per hour brackets the empirical breakpoints found by
#'   Faber et al. (2023) and Jørgensen et al. (2022) across taxa as different
#'   as *Drosophila suzukii* and *Lemna gibba*.
#'
#' All three quantities inherit the same posterior — no extra fitting step.
#'
#' @param workflow    Fitted `tdt_4pl_workflow`.
#' @param t_ref       Reference exposure duration for CTmax, in the
#'                    `output_time_unit` (default `"min"`). Default 60.
#' @param TC_rate_range Numeric length-2: HI-rate floor range, in % LT50-dose
#'                    per hour, used to derive T_crit. Default `c(0.1, 1)`.
#'                    Sampled uniformly on `log10(r/100)`, which is the natural
#'                    scale for a rate threshold.
#' @param temp_grid   Numeric vector of temperatures to search over. Default:
#'                    a fine grid spanning the training-data temperature range
#'                    extended by ±2 °C.
#' @param duration_grid Optional duration grid for the underlying LT50 curve.
#'                    Default: 350 log-spaced values spanning 0.2× to 5× the
#'                    training-data duration range.
#' @param ndraws      Posterior draws to use. Default 1000.
#' @param time_multiplier Multiplier from model time units to `output_time_unit`.
#'                    Default 60 (so hour-scale data → minute-scale outputs).
#' @param output_time_unit Label for the output time unit. Default `"min"`.
#' @return A list with elements:
#'   - `z`: list with `draws` (per-draw z, slope, intercept, R²) and `summary`.
#'   - `CTmax`: list with `draws` (per-draw temperature) and `summary`.
#'   - `T_crit`: list with `draws` and `summary`.
#'   - `lt50_curve`: output of [derive_ltx_curve()] (intermediate).
#'   - `meta`: list of inputs used (`t_ref`, `TC_rate_range`, `output_time_unit`).
#' @examples
#' \dontrun{
#' wf  <- fit_4pl(d, ...)
#' out <- extract_tdt(wf)
#' out$z$summary
#' out$CTmax$summary
#' out$T_crit$summary
#' # Feed the T_crit posterior median into predict_heat_injury():
#' hi <- predict_heat_injury(trace, wf, T_c = out$T_crit$summary$temp_median)
#' }
#' @export
extract_tdt <- function(workflow,
                        t_ref            = 60,
                        TC_rate_range    = c(0.1, 1),
                        temp_grid        = NULL,
                        duration_grid    = NULL,
                        ndraws           = 1000,
                        time_multiplier  = 60,
                        output_time_unit = "min") {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  if (length(TC_rate_range) != 2L ||
      any(!is.finite(TC_rate_range)) ||
      any(TC_rate_range <= 0) ||
      TC_rate_range[1] >= TC_rate_range[2])
    stop("TC_rate_range must be c(low, high) with 0 < low < high (% HI/hour).",
         call. = FALSE)

  data <- workflow$data
  if (is.null(temp_grid)) {
    trange   <- range(data$temp, na.rm = TRUE)
    temp_grid <- seq(trange[1] - 2, trange[2] + 2, by = 0.05)
  }

  # LT50 curve → z + CTmax via per-draw log-linear regression
  lt50_curve <- derive_ltx_curve(
    workflow         = workflow,
    temp_grid        = temp_grid,
    duration_grid    = duration_grid,
    target_surv      = 0.5,
    ndraws           = ndraws,
    time_multiplier  = time_multiplier,
    output_time_unit = output_time_unit
  )
  z_ctmax <- derive_tdt_parameters(lt50_curve, t_ref = t_ref)

  # Express the t_ref reference duration back in model time units so the
  # inverse-4PL lookup uses the same scale as the fitted model.
  exposure_in_model_units <- t_ref / time_multiplier

  ctmax <- derive_temperature_for_duration(
    workflow          = workflow,
    exposure_duration = exposure_in_model_units,
    temp_grid         = temp_grid,
    target_surv       = 0.5,
    ndraws            = ndraws
  )

  # T_crit via rate-multiplier integration. For each posterior draw, sample
  # r* uniformly on log10 across TC_rate_range, then compute
  # T_crit = CTmax + z * log10(r*/100). Pairs z and CTmax by .draw index
  # so the joint posterior is preserved; inner_join drops draws that didn't
  # survive both pipelines (e.g. ill-conditioned LT50 regressions).
  z_df     <- z_ctmax$draws   |> dplyr::select(.draw, z)
  ctmax_df <- ctmax$draws     |> dplyr::select(.draw, CTmax_temp = temp)
  paired   <- dplyr::inner_join(z_df, ctmax_df, by = ".draw")

  log10_low  <- log10(TC_rate_range[1] / 100)
  log10_high <- log10(TC_rate_range[2] / 100)
  paired$log10_rate <- stats::runif(nrow(paired),
                                     min = log10_low, max = log10_high)
  paired$T_crit     <- paired$CTmax_temp + paired$z * paired$log10_rate

  t_crit_draws <- tibble::tibble(.draw = paired$.draw,
                                  temp = paired$T_crit,
                                  log10_rate = paired$log10_rate)
  q <- stats::quantile(t_crit_draws$temp,
                       c(0.025, 0.5, 0.975), na.rm = TRUE, names = FALSE)
  t_crit_summary <- tibble::tibble(
    TC_rate_low  = TC_rate_range[1],
    TC_rate_high = TC_rate_range[2],
    temp_lower   = q[1],
    temp_median  = q[2],
    temp_upper   = q[3]
  )

  list(
    z          = list(draws   = z_ctmax$draws |> dplyr::select(.draw, z),
                      summary = z_ctmax$summary |>
                        dplyr::select(z_median, z_lower, z_upper)),
    CTmax      = list(draws   = ctmax$draws,
                      summary = ctmax$summary),
    T_crit     = list(draws   = t_crit_draws,
                      summary = t_crit_summary),
    lt50_curve = lt50_curve,
    meta       = list(t_ref            = t_ref,
                      TC_rate_range    = TC_rate_range,
                      output_time_unit = output_time_unit)
  )
}
