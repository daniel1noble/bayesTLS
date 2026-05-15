# Classical TDT quantities (z, CTmax, T_crit) derived from a fitted 4PL via
# either the 4PL `mid` parameter directly (relative threshold; default) or
# numerical inversion of the posterior survival surface at an absolute
# survival probability.

#' Normalise a `target_surv` argument
#'
#' Accepts the user-facing argument (string `"relative"`/`"absolute"` or
#' numeric in `(0, 1)`) and returns a list describing the chosen threshold
#' mode plus a character label suitable for embedding in result tibbles.
#'
#' - `"relative"` (default) → threshold = `(low + up)/2` per posterior draw.
#'   The 4PL `mid` parameter is the log10-time at this threshold, so no
#'   numerical inversion is needed.
#' - `"absolute"` → threshold = 0.5 (literal survival probability).
#' - numeric `p` in `(0, 1)` → threshold = `p` (literal survival probability).
#'
#' @keywords internal
resolve_target_surv <- function(target_surv) {
  if (is.character(target_surv) && length(target_surv) == 1L) {
    if (target_surv == "relative") {
      return(list(mode = "relative", prob = NA_real_,
                  label = "(low+up)/2"))
    }
    if (target_surv == "absolute") {
      return(list(mode = "absolute", prob = 0.5,
                  label = sprintf("p=%.3f", 0.5)))
    }
    stop("target_surv must be \"relative\", \"absolute\", or a numeric in (0, 1).",
         call. = FALSE)
  }
  if (is.numeric(target_surv) && length(target_surv) == 1L &&
      is.finite(target_surv) && target_surv > 0 && target_surv < 1) {
    return(list(mode = "absolute", prob = as.numeric(target_surv),
                label = sprintf("p=%.3f", as.numeric(target_surv))))
  }
  stop("target_surv must be \"relative\", \"absolute\", or a numeric in (0, 1).",
       call. = FALSE)
}

#' Posterior LT_x curve: time to reach a survival target at each temperature
#'
#' Returns the per-draw duration at which population-level survival crosses
#' the chosen threshold, at each temperature in `temp_grid`.
#'
#' Two threshold modes are supported via `target_surv`:
#'
#' - `"relative"` (default): the duration at which survival reaches the
#'   midpoint between the fitted lower and upper asymptotes, i.e.
#'   `(low + up)/2`. This is the 4PL `mid` parameter on the natural time
#'   axis, returned directly from `posterior_linpred(nlpar = "mid")` — no
#'   numerical inversion. When `low ≈ 0` and `up ≈ 1` it coincides with
#'   the classical LT50.
#' - `"absolute"` (or a numeric `p` in `(0, 1)`): the duration at which
#'   survival crosses the literal probability `p` (0.5 by default). The
#'   inversion is numerical — predict survival on a dense duration grid,
#'   then `approx()` through `p`.
#'
#' This is the **horizontal** read of the survival surface: fix a survival
#' threshold, read off the time required to reach it at each temperature.
#'
#' @param workflow         Fitted `bayes_tls`.
#' @param temp_grid        Numeric vector of temperatures (°C).
#' @param duration_grid    Numeric vector of durations along which to search.
#'                         Only used in `"absolute"` mode. Default: 350
#'                         log-spaced values spanning 0.2× to 5× the training
#'                         data's duration range.
#' @param target_surv      Threshold mode. `"relative"` (default), `"absolute"`
#'                         (= 0.5), or a numeric in `(0, 1)`.
#' @param ndraws           Posterior draws to use. Default 1000.
#' @param probs            Quantile probabilities for the summary. Default
#'                         `c(0.025, 0.5, 0.975)`.
#' @param time_multiplier  Multiplier from model time units to output time
#'                         units (e.g. 60 for hours → min). Default 60.
#' @param output_time_unit Label for the output time unit. Default `"min"`.
#' @return A list with `draws` (per-draw threshold durations; `target_surv`
#'         column is a character label), `summary` (quantile summary by
#'         temperature), `target_surv` (the label), `time_multiplier`,
#'         `output_time_unit`.
#' @export
derive_tdt_curve <- function(workflow,
                             temp_grid,
                             duration_grid    = NULL,
                             target_surv      = "relative",
                             ndraws           = 1000,
                             probs            = c(0.025, 0.5, 0.975),
                             time_multiplier  = 60,
                             output_time_unit = "min") {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  ts <- resolve_target_surv(target_surv)

  if (ts$mode == "relative") {
    # Direct shortcut: log10(t_relative) = mid(T) per draw. No grid search.
    nd <- new_tdt_grid(workflow, temps = temp_grid, durations = 1)
    pp_mid <- brms::posterior_linpred(workflow$fit, newdata = nd,
                                       nlpar = "mid", re_formula = NA,
                                       ndraws = ndraws)
    # pp_mid is [ndraws x length(temp_grid)] of log10(t) in model time units.
    duration_model_mat <- 10 ^ pp_mid

    draw_list <- vector("list", length(temp_grid))
    for (i in seq_along(temp_grid)) {
      t_i  <- temp_grid[i]
      dmod <- duration_model_mat[, i]
      draw_list[[i]] <- data.frame(
        .draw            = seq_along(dmod),
        temp             = t_i,
        target_surv      = ts$label,
        duration_model   = dmod,
        duration_out     = dmod * time_multiplier,
        stringsAsFactors = FALSE
      )
    }
  } else {
    if (is.null(duration_grid)) {
      drange <- range(workflow$data$duration, na.rm = TRUE)
      duration_grid <- 10 ^ seq(log10(drange[1] / 5),
                                log10(drange[2] * 5),
                                length.out = 350)
    }
    nd   <- new_tdt_grid(workflow, temps = temp_grid, durations = duration_grid)
    pred <- posterior_linpred_tdt(workflow, nd, ndraws = ndraws,
                                   re_formula = NA)
    draw_list <- vector("list", length(temp_grid))
    for (i in seq_along(temp_grid)) {
      t_i <- temp_grid[i]
      idx <- nd$temp == t_i
      thr <- threshold_x_by_draw(pred_mat = pred[, idx, drop = FALSE],
                                 x        = nd$duration[idx],
                                 target   = ts$prob)
      draw_list[[i]] <- data.frame(
        .draw            = seq_along(thr),
        temp             = t_i,
        target_surv      = ts$label,
        duration_model   = thr,
        duration_out     = thr * time_multiplier,
        stringsAsFactors = FALSE
      )
    }
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
       target_surv      = ts$label,
       target_mode      = ts$mode,
       target_prob      = ts$prob,
       time_multiplier  = time_multiplier,
       output_time_unit = output_time_unit)
}

#' Temperature at which survival equals a target after a fixed exposure
#'
#' The **vertical** read of the survival surface: fix an exposure duration,
#' find the temperature at which the posterior survival reaches the chosen
#' threshold. Returns one temperature per posterior draw.
#'
#' Threshold modes (via `target_surv`) match [derive_tdt_curve()]:
#'
#' - `"relative"` (default) → temperature at which `mid(T) = log10(exposure_duration)`
#'   per draw. The inversion is done analytically per draw: extract
#'   `posterior_linpred(nlpar = "mid")` over `temp_grid`, then `approx()` to
#'   the target log10-time.
#' - `"absolute"` (= 0.5) or numeric `p` in `(0, 1)` → existing numerical
#'   inversion of the 4PL survival surface at the literal probability `p`.
#'
#' This is the primitive used by [extract_tdt()] to derive CTmax at `t_ref`.
#'
#' @param workflow         Fitted `bayes_tls`.
#' @param exposure_duration Numeric scalar — the fixed duration (model units).
#' @param temp_grid        Numeric vector of temperatures to search over.
#' @param target_surv      Threshold mode. `"relative"` (default), `"absolute"`,
#'                         or a numeric in `(0, 1)`.
#' @param ndraws           Posterior draws. Default 1000.
#' @param probs            Quantile probabilities. Default `c(0.025, 0.5, 0.975)`.
#' @return A list with `draws` (per-draw threshold temperatures; `target_surv`
#'         column is a character label), `summary` (quantile summary),
#'         `exposure_duration`, `target_surv` (the label).
#' @export
derive_temperature_for_duration <- function(workflow,
                                            exposure_duration,
                                            temp_grid,
                                            target_surv = "relative",
                                            ndraws      = 1000,
                                            probs       = c(0.025, 0.5, 0.975)) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  ts <- resolve_target_surv(target_surv)

  if (ts$mode == "relative") {
    # Per-draw inversion of mid(T) = log10(exposure_duration) in model time.
    nd_mid <- new_tdt_grid(workflow, temps = temp_grid, durations = 1)
    pp_mid <- brms::posterior_linpred(workflow$fit, newdata = nd_mid,
                                       nlpar = "mid", re_formula = NA,
                                       ndraws = ndraws)
    target_logd <- log10(exposure_duration)
    # pp_mid is [ndraws x length(temp_grid)] of mid(T) per draw.
    thr <- vapply(seq_len(nrow(pp_mid)), function(d) {
      mid_vec <- pp_mid[d, ]
      if (!any(is.finite(mid_vec))) return(NA_real_)
      ord <- order(mid_vec)
      suppressWarnings(
        stats::approx(mid_vec[ord], temp_grid[ord], xout = target_logd)$y
      )
    }, numeric(1))
  } else {
    nd   <- new_tdt_grid(workflow, temps = temp_grid,
                         durations = exposure_duration)
    pred <- posterior_linpred_tdt(workflow, nd, ndraws = ndraws,
                                   re_formula = NA)
    thr  <- threshold_x_by_draw(pred_mat = pred,
                                x        = nd$temp,
                                target   = ts$prob)
  }

  draws <- tibble::tibble(
    .draw       = seq_along(thr),
    target_surv = ts$label,
    temp        = thr
  ) |>
    dplyr::filter(is.finite(temp))

  summary <- draws |>
    dplyr::summarise(
      target_surv = ts$label,
      temp_lower  = stats::quantile(temp, probs[1], na.rm = TRUE),
      temp_median = stats::quantile(temp, probs[2], na.rm = TRUE),
      temp_upper  = stats::quantile(temp, probs[3], na.rm = TRUE)
    )

  list(draws             = draws,
       summary           = summary,
       exposure_duration = exposure_duration,
       target_surv       = ts$label,
       target_mode       = ts$mode,
       target_prob       = ts$prob)
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
#' @param ltx_curve  Output of [derive_tdt_curve()].
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

#' Extract classical TDT quantities from a fitted 4PL: z, CTmax (optional T_crit)
#'
#' Always returns:
#'
#' - **z** — thermal sensitivity, derived by per-draw log-linear regression of
#'   `log10(LT_x)` on temperature (see [derive_tdt_parameters()]). The
#'   threshold defining the LT_x curve is controlled by `target_surv`.
#' - **CTmax** — temperature at which survival reaches the chosen threshold
#'   after `t_ref` exposure; the temperature where the LT_x curve crosses
#'   `t_ref`.
#'
#' By default (`target_surv = "relative"`), the threshold is the per-draw
#' midpoint between the fitted lower and upper asymptotes (`(low + up)/2`).
#' This is the most biologically meaningful threshold when the upper asymptote
#' is below 1 (e.g., when there is intrinsic background mortality unrelated to
#' heat stress). When `low ≈ 0` and `up ≈ 1`, this coincides with the
#' classical absolute 50 % LT50. Pass `target_surv = "absolute"` (or any
#' numeric in `(0, 1)`) to recover the absolute threshold.
#'
#' When `lethal = TRUE` it *also* returns **T_crit**, the rate-multiplier
#' critical temperature: for each posterior draw,
#' `T_crit = CTmax + z * log10(r* / 100)`, with `r*` drawn uniformly on the
#' `log10` scale across `TC_rate_range`. The pooled posterior thus carries
#' both parameter uncertainty (in `CTmax` and `z`) and operational uncertainty
#' in the choice of damage-rate floor. The default range `c(0.1, 1)` %
#' HI per hour brackets the empirical breakpoints found by Faber et al. (2026)
#' and Jørgensen et al. (2021) across taxa as different as *Drosophila suzukii*
#' and *Lemna gibba*.
#'
#' T_crit only makes physical sense for **lethal-endpoint** data — proportion-
#' or count-based survival under a damage-accumulation interpretation. For
#' sublethal endpoints (knockdown time, photosystem-II failure, etc.) the
#' fitted `z` measures the rate of *performance reduction* rather than damage
#' accumulation, and the two are not interchangeable: sublethal `z` is
#' typically far steeper, which in turn pushes the rate-multiplier `T_crit`
#' implausibly low. Setting `lethal = FALSE` (the default) suppresses `T_crit`
#' to avoid that pitfall; users with lethal data opt in by passing
#' `lethal = TRUE` and are reminded by a startup message.
#'
#' @param workflow    Fitted `bayes_tls`.
#' @param target_surv Threshold mode. `"relative"` (default; threshold =
#'                    `(low + up)/2` per draw, computed from `mid`),
#'                    `"absolute"` (= 0.5), or a numeric in `(0, 1)`.
#' @param t_ref       Reference exposure duration for CTmax, in the
#'                    `output_time_unit` (default `"min"`). Default 60.
#' @param TC_rate_range Numeric length-2: HI-rate floor range, in % LT-dose
#'                    per hour, used to derive T_crit (only when
#'                    `lethal = TRUE`). Default `c(0.1, 1)`. Sampled uniformly
#'                    on `log10(r/100)`, which is the natural scale for a
#'                    rate threshold.
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
#' @param lethal      Logical. When `TRUE`, also returns the rate-multiplier
#'                    T_crit and emits a one-line reminder that T_crit is
#'                    valid only for damage-accumulation (lethal) endpoints.
#'                    Default `FALSE`.
#' @return A list with elements:
#'   - `z`: list with `draws` (per-draw z, slope, intercept, R²) and `summary`.
#'   - `CTmax`: list with `draws` (per-draw temperature) and `summary`.
#'   - `T_crit`: list with `draws` and `summary` when `lethal = TRUE`;
#'     `NULL` otherwise.
#'   - `lt50_curve`: output of [derive_tdt_curve()] (intermediate).
#'   - `meta`: list of inputs used (`t_ref`, `TC_rate_range`, `lethal`,
#'     `output_time_unit`).
#' @examples
#' \dontrun{
#' wf  <- fit_4pl(d, ...)
#' out <- extract_tdt(wf)                  # z + CTmax only
#' out$z$summary
#' out$CTmax$summary
#'
#' # Lethal-endpoint data — opt in to T_crit:
#' out2 <- extract_tdt(wf, lethal = TRUE)
#' out2$T_crit$summary
#' # Feed the T_crit posterior median into predict_heat_injury():
#' hi <- predict_heat_injury(trace, wf, T_c = out2$T_crit$summary$temp_median)
#' }
#' @export
extract_tdt <- function(workflow,
                        target_surv      = "relative",
                        t_ref            = 60,
                        TC_rate_range    = c(0.1, 1),
                        temp_grid        = NULL,
                        duration_grid    = NULL,
                        ndraws           = 1000,
                        time_multiplier  = 60,
                        output_time_unit = "min",
                        lethal           = FALSE) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  if (length(TC_rate_range) != 2L ||
      any(!is.finite(TC_rate_range)) ||
      any(TC_rate_range <= 0) ||
      TC_rate_range[1] >= TC_rate_range[2])
    stop("TC_rate_range must be c(low, high) with 0 < low < high (% HI/hour).",
         call. = FALSE)

  # Validate up front so both helpers receive a normalised label.
  ts <- resolve_target_surv(target_surv)

  data <- workflow$data
  if (is.null(temp_grid)) {
    trange   <- range(data$temp, na.rm = TRUE)
    temp_grid <- seq(trange[1] - 2, trange[2] + 2, by = 0.05)
  }

  # LT_x curve → z + CTmax via per-draw log-linear regression. The threshold
  # for the curve is set by target_surv; `"relative"` (default) returns mid(T)
  # directly.
  lt50_curve <- derive_tdt_curve(
    workflow         = workflow,
    temp_grid        = temp_grid,
    duration_grid    = duration_grid,
    target_surv      = target_surv,
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
    target_surv       = target_surv,
    ndraws            = ndraws
  )

  t_crit_block <- NULL
  if (isTRUE(lethal)) {
    message("extract_tdt(): T_crit reported under the rate-multiplier ",
            "definition; valid for damage-accumulation (lethal) endpoints ",
            "only. If your data are sublethal (knockdown, performance, ",
            "PSII, ...) ignore T_crit and supply T_c manually downstream.")

    # T_crit via rate-multiplier integration. For each posterior draw, sample
    # r* uniformly on log10 across TC_rate_range, then compute
    # T_crit = CTmax + z * log10(r*/100). Pairs z and CTmax by .draw index
    # so the joint posterior is preserved.
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
    t_crit_block <- list(draws = t_crit_draws, summary = t_crit_summary)
  }

  list(
    z          = list(draws   = z_ctmax$draws |> dplyr::select(.draw, z),
                      summary = z_ctmax$summary |>
                        dplyr::select(z_median, z_lower, z_upper)),
    CTmax      = list(draws   = ctmax$draws,
                      summary = ctmax$summary),
    T_crit     = t_crit_block,
    lt50_curve = lt50_curve,
    meta       = list(target_surv      = ts$label,
                      target_mode      = ts$mode,
                      target_prob      = ts$prob,
                      t_ref            = t_ref,
                      TC_rate_range    = TC_rate_range,
                      lethal           = isTRUE(lethal),
                      output_time_unit = output_time_unit)
  )
}
