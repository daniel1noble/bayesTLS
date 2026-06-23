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
#' @param target_surv      Threshold mode. `"relative"` (default; = `(low + up)/2`),
#'                         `"absolute"` (= 0.5), or a numeric in `(0, 1)`.
#' @param ndraws           Posterior draws to use. Default 1000.
#' @param probs            Quantile probabilities for the summary. Default
#'                         `c(0.025, 0.5, 0.975)`.
#' @param time_multiplier  Multiplier from model time units to `output_time_unit`
#'                         (e.g. 60 for an hours model → min). `NULL` (default)
#'                         derives it automatically from the workflow's
#'                         `duration_unit` and `output_time_unit`, so a minutes
#'                         model and an hours model both give the correct result
#'                         without manual tuning. Pass a value to override.
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
                             time_multiplier  = NULL,
                             output_time_unit = "min") {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  # brms::posterior_linpred() errors if `ndraws` exceeds the posterior size, so
  # clamp here (the default ndraws = 1000 otherwise crashes the relative-mode
  # call below on any fit with fewer draws, e.g. 2 chains x 400). The absolute
  # branch goes through posterior_linpred_tdt(), which clamps the same way.
  if (!is.null(ndraws)) {
    total <- tryCatch(brms::ndraws(workflow$fit), error = function(e) NA_integer_)
    if (is.finite(total)) ndraws <- min(ndraws, total)
  }

  time_multiplier <- tdt_resolve_time_multiplier(time_multiplier, workflow$meta,
                                                 output_time_unit)
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
#' @param ndraws           Posterior draws to subsample, or `NULL` for all.
#'                         Default 1000.
#' @param probs            Quantile probabilities. Default `c(0.025, 0.5, 0.975)`.
#' @param seed             Optional integer seeding the draw subsample for
#'                         reproducibility. `NULL` (default) leaves the RNG alone.
#' @param by               Optional moderator column(s) for per-group results.
#'                         `NULL` (default) uses the fit's moderators; a
#'                         single-condition fit returns one ungrouped result.
#' @return A list with `draws` (per-draw threshold temperatures; `target_surv`
#'         column is a character label), `summary` (quantile summary),
#'         `exposure_duration`, `target_surv` (the label), `target_mode`,
#'         `target_prob`. A grouped fit adds the moderator column(s).
#' @export
derive_temperature_for_duration <- function(workflow,
                                            exposure_duration,
                                            temp_grid,
                                            target_surv = "absolute",
                                            ndraws      = 1000,
                                            probs       = c(0.025, 0.5, 0.975),
                                            seed        = NULL,
                                            by          = NULL) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  ts     <- resolve_target_surv(target_surv)
  by     <- tdt_resolve_by(workflow, by)
  tbar   <- workflow$meta$temp_mean
  target <- log10(exposure_duration)
  if (!is.null(seed)) set.seed(seed)   # reproducible draw subsample
  fit <- get_brmsfit(workflow)
  did <- tls_draw_ids(fit, ndraws)

  # Relative: closed-form inverse of the linear midpoint (2-temp grid). Absolute:
  # numerical inversion of the bent logLT curve over the search grid.
  if (ts$mode == "relative") {
    nd <- tls_build_grid(fit$data, by = by, temp = "temp_c", temp_grid = c(0, 1))
    sp <- tls_eval_subpars(fit, nd, workflow$meta$bounds, draw_ids = did, mode = "relative")
  } else {
    nd <- tls_build_grid(fit$data, by = by, temp = "temp_c", temp_grid = temp_grid - tbar)
    sp <- tls_eval_subpars(fit, nd, workflow$meta$bounds, draw_ids = did,
                           mode = "absolute", p = ts$prob %||% 0.5)
  }

  per_group <- lapply(unique(nd$.grp), function(g) {
    gi <- which(nd$.grp == g)
    Tc <- if (ts$mode == "relative") {
      c0 <- gi[which(nd$temp_c[gi] == 0)]; c1 <- gi[which(nd$temp_c[gi] == 1)]
      tbar + (target - sp$mid[, c0]) / (sp$mid[, c1] - sp$mid[, c0])
    } else {
      tls_invert_logLT(sp$logLT[, gi, drop = FALSE], target, temp_grid)
    }
    d <- tibble::tibble(.draw = seq_along(Tc), target_surv = ts$label, temp = Tc) |>
      dplyr::filter(is.finite(temp))
    q <- stats::quantile(d$temp, probs, names = FALSE, na.rm = TRUE)
    s <- tibble::tibble(target_surv = ts$label, temp_lower = q[1],
                        temp_median = q[2], temp_upper = q[3])
    if (!is.null(by)) {
      gc <- nd[gi[1], by, drop = FALSE]
      d <- cbind(gc, d, row.names = NULL); s <- cbind(gc, s, row.names = NULL)
    }
    list(draws = d, summary = s)
  })

  list(draws             = dplyr::bind_rows(lapply(per_group, `[[`, "draws")),
       summary           = dplyr::bind_rows(lapply(per_group, `[[`, "summary")),
       exposure_duration = exposure_duration,
       target_surv       = ts$label,
       target_mode       = ts$mode,
       target_prob       = ts$prob)
}

#' Per-draw thermal sensitivity z directly from the joint posterior
#'
#' Derives \eqn{z = -1 / (\mathrm{d}/\mathrm{d}T\,\log_{10}\mathrm{LT}(T))} per
#' posterior draw, read straight from the fitted 4PL coefficients — **no
#' regression**. There are two regimes:
#'
#' - **Relative threshold** (default; the \eqn{(\ell+u)/2} midpoint):
#'   \eqn{\log_{10}\mathrm{LT}_{\text{rel}}(T) = \mathrm{mid}(T) =
#'   \beta_0 + \beta_1 (T-\bar T)} is exactly linear, so
#'   \eqn{z = -1/\beta_1} where \eqn{\beta_1} is the temperature slope on
#'   `mid` (`b_mid_temp_c`). The asymptotes \eqn{\ell, u} and slope \eqn{k}
#'   do not enter (the midpoint cancels the curve asymmetry). z is constant in
#'   temperature.
#' - **Absolute threshold** \eqn{p}: the LT curve gains the asymmetry-correction
#'   term, \eqn{\log_{10}\mathrm{LT}_p(T) = \mathrm{mid}(T) +
#'   \tfrac{1}{k(T)}\log\tfrac{u(T)-p}{p-\ell(T)}}. When \eqn{\ell}, \eqn{u} or
#'   \eqn{k} carry temperature effects this bends the curve, so z varies with
#'   temperature. A **local** \eqn{z(T) = -1/m(T)} is computed at each assay
#'   temperature, where the local slope \eqn{m(T)} is obtained by a central
#'   finite difference of the closed-form LT curve (step `h`). When the shape
#'   parameters are constant in T the correction is flat and this reduces to
#'   \eqn{-1/\beta_1}.
#'
#' The returned **pooled** z (the default single summary) is the per-draw mean
#' of the local \eqn{z(T)} over `temp_grid`. The full per-temperature local
#' \eqn{z(T)} is also returned. z is invariant to the time unit (a constant
#' time-multiplier shifts the LT intercept, not its slope), so no
#' `time_multiplier` is needed here.
#'
#' @param workflow    Fitted `bayes_tls`.
#' @param target_surv Threshold mode: `"relative"` (default; = `(low + up)/2`),
#'                    `"absolute"` (= 0.5), or a numeric in `(0, 1)`.
#' @param temp_grid   Temperatures at which to evaluate local z and over which
#'                    to pool. Default: the observed (unique) assay temperatures
#'                    — pooling only where the data inform the curve.
#' @param ndraws      Posterior draws to subsample, or `NULL` (default) for all.
#' @param probs       Quantile probabilities for the summaries. Default
#'                    `c(0.025, 0.5, 0.975)`.
#' @param h           Temperature step (°C) for the central finite difference.
#'                    Default `1e-3`. (For a linear midpoint — the relative
#'                    threshold — the central difference is exact regardless.)
#' @param seed        Optional integer seeding the draw subsample (relevant only
#'                    when `ndraws` is set) for reproducibility. `NULL` (default)
#'                    leaves the RNG untouched.
#' @param by          Optional moderator column(s) for per-group z. `NULL`
#'                    (default) uses the fit's moderators (`meta$group_vars`); a
#'                    single-condition fit then returns one ungrouped result.
#' @return A list with:
#'   - `draws`: tibble `(.draw, z)` — pooled per-draw z.
#'   - `summary`: tibble `(z_median, z_lower, z_upper)`.
#'   - `local_draws`: tibble `(.draw, temp, z)` — local z(T) per draw.
#'   - `local_summary`: tibble `(temp, z_median, z_lower, z_upper)`.
#'   - `target_surv`, `temp_grid`.
#'   For a grouped fit each tibble gains the moderator column(s).
#' @examples
#' \dontrun{
#' wf <- fit_4pl(std)
#' z  <- derive_z(wf)             # relative: z = -1 / b_mid_temp_c per draw
#' z$summary
#' derive_z(wf, target_surv = "absolute")$local_summary  # local z(T)
#' }
#' @export
derive_z <- function(workflow,
                     target_surv = "relative",
                     temp_grid   = NULL,
                     ndraws      = NULL,
                     probs       = c(0.025, 0.5, 0.975),
                     h           = 1e-3,
                     seed        = NULL,
                     by          = NULL) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  ts   <- resolve_target_surv(target_surv)
  by   <- tdt_resolve_by(workflow, by)
  tbar <- workflow$meta$temp_mean

  if (is.null(temp_grid)) temp_grid <- sort(unique(workflow$data$temp))
  temp_grid <- temp_grid[is.finite(temp_grid)]
  if (length(temp_grid) < 1L) stop("temp_grid is empty.", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)   # reproducible draw subsample (when ndraws set)
  fit  <- get_brmsfit(workflow)
  did  <- tls_draw_ids(fit, ndraws)
  L    <- length(temp_grid)
  tc   <- temp_grid - tbar
  # Evaluate logLT at temp_grid - h (first L cols per group) and + h (next L), so
  # z(T) = -1 / central-difference. Linear (relative) mid -> exact; absolute ->
  # local slope of the bent curve. Same maths as before, posterior_linpred source.
  nd <- tls_build_grid(fit$data, by = by, temp = "temp_c", temp_grid = c(tc - h, tc + h))
  sp <- tls_eval_subpars(fit, nd, workflow$meta$bounds, draw_ids = did,
                         mode = if (ts$mode == "relative") "relative" else "absolute",
                         p = ts$prob %||% 0.5)

  per_group <- lapply(unique(nd$.grp), function(g) {
    gi <- which(nd$.grp == g)                      # 2L cols: 1:L minus, (L+1):2L plus
    zo <- tls_local_z(sp$logLT[, gi[(L + 1):(2 * L)], drop = FALSE],
                      sp$logLT[, gi[1:L], drop = FALSE], h, temp_grid, probs)
    if (!is.null(by)) {
      gc <- nd[gi[1], by, drop = FALSE]
      zo$draws         <- cbind(gc, zo$draws,         row.names = NULL)
      zo$summary       <- cbind(gc, zo$summary,       row.names = NULL)
      zo$local_draws   <- cbind(gc, zo$local_draws,   row.names = NULL)
      zo$local_summary <- cbind(gc, zo$local_summary, row.names = NULL)
    }
    zo
  })
  cmb <- if (is.null(by)) per_group[[1]] else list(
    draws         = dplyr::bind_rows(lapply(per_group, `[[`, "draws")),
    summary       = dplyr::bind_rows(lapply(per_group, `[[`, "summary")),
    local_draws   = dplyr::bind_rows(lapply(per_group, `[[`, "local_draws")),
    local_summary = dplyr::bind_rows(lapply(per_group, `[[`, "local_summary")))
  c(cmb, list(target_surv = ts$label, temp_grid = temp_grid))
}

#' Extract classical TDT quantities from a fitted 4PL: z, CTmax (optional T_crit)
#'
#' Always returns:
#'
#' - **z** — thermal sensitivity, read directly from the joint posterior (no
#'   regression): the relative threshold gives `z = -1 / b_mid_temp_c` per
#'   draw; an absolute threshold pools the per-draw local `z(T)` over the assay
#'   temperatures (see [derive_z()]). The threshold is controlled by
#'   `target_surv`.
#' - **CTmax** — temperature at which survival reaches the chosen threshold
#'   after `t_ref` exposure; the temperature where the LT_x curve crosses
#'   `t_ref`, by per-draw inversion of the fitted surface.
#'
#' By default (`target_surv = "relative"`), the threshold is the per-draw
#' midpoint between the fitted lower and upper asymptotes (`(low + up)/2`); it
#' is more biologically meaningful when the upper asymptote is below 1 (e.g.,
#' intrinsic background mortality unrelated to heat stress), and coincides with
#' the absolute 50 % LT50 when `low ≈ 0` and `up ≈ 1`. Pass
#' `target_surv = "absolute"` for the field-standard absolute 50 % survival
#' level (LT50), or any numeric in `(0, 1)` for a custom threshold.
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
#'                    `(low + up)/2` per draw, from `mid`), `"absolute"`
#'                    (the 50 % LT50), or a numeric in `(0, 1)`.
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
#'                    `NULL` (default) derives it from the workflow's
#'                    `duration_unit` and `output_time_unit` (e.g. 60 for an
#'                    hours model, 1 for a minutes model). Pass a value to
#'                    override. This is what makes `t_ref` (in `output_time_unit`)
#'                    map to the correct exposure regardless of the model's time
#'                    unit — omitting it on a minutes model used to compute CTmax
#'                    at `t_ref/60`.
#' @param output_time_unit Label for the output time unit. Default `"min"`.
#' @param lethal      Logical. When `TRUE`, also returns the rate-multiplier
#'                    T_crit and emits a one-line reminder that T_crit is
#'                    valid only for damage-accumulation (lethal) endpoints.
#'                    Default `FALSE`.
#' @param z_local     Logical. When `TRUE`, additionally computes and returns
#'                    the per-draw local `z(T)` at each assay temperature in
#'                    `z$local` (relevant when an absolute threshold and
#'                    temperature-varying asymptotes bend the LT curve). When
#'                    `FALSE` (default) this per-temperature breakdown is not
#'                    computed, which saves the dominant per-call cost; the
#'                    pooled `z` draws and summary are identical either way.
#'                    For a model R\eqn{^2}, call
#'                    `brms::bayes_R2(get_brmsfit(workflow))`.
#' @param seed        Optional integer. When supplied, seeds the RNG so the draw
#'                    subsample and the `T_crit` rate draws are reproducible.
#'                    `NULL` (default) leaves the RNG untouched.
#' @param by          Optional moderator column(s) for per-group results. `NULL`
#'                    (default) uses the fit's moderators (`meta$group_vars`); a
#'                    single-condition fit then returns the ungrouped result
#'                    (the original nested shape). For a grouped fit, `z`,
#'                    `CTmax` and `T_crit` summaries/draws gain the moderator
#'                    column(s) (one block per group); `lt50_curve` is `NULL`.
#' @return A list with elements:
#'   - `z`: list with `draws` (per-draw pooled z) and `summary`; plus `local`
#'     (`draws` + per-temperature `summary`) when `z_local = TRUE`, else `NULL`.
#'     z is read directly from the posterior — relative threshold gives
#'     `-1 / b_mid_temp_c` per draw; an absolute threshold pools the per-draw
#'     local `z(T)` over the assay temperatures (see [derive_z()]).
#'   - `CTmax`: list with `draws` (per-draw temperature) and `summary`.
#'   - `T_crit`: list with `draws` and `summary` when `lethal = TRUE`;
#'     `NULL` otherwise. z and CTmax share the same posterior draws, so the
#'     pairing is genuinely joint.
#'   - `lt50_curve`: output of [derive_tdt_curve()] (descriptive intermediate).
#'   - `meta`: list of inputs used (`t_ref`, `TC_rate_range`, `lethal`,
#'     `output_time_unit`).
#' @seealso [get_z_summary()], [get_ctmax_summary()], [get_tcrit_summary()] for
#'   tidy summary tibbles, and [get_z_draws()], [get_ctmax_draws()],
#'   [get_tcrit_draws()] for the per-draw posteriors (e.g. group contrasts).
#' @examples
#' \dontrun{
#' wf  <- fit_4pl(d, ...)
#' out <- extract_tdt(wf)                  # z + CTmax only
#' get_z_summary(out)                      # median + 95% CrI
#' get_ctmax_summary(out)
#' get_z_draws(out)                        # per-draw posterior (e.g. contrasts)
#'
#' # Lethal-endpoint data — opt in to T_crit:
#' out2 <- extract_tdt(wf, lethal = TRUE)
#' get_tcrit_summary(out2)
#' # Feed the T_crit posterior median into predict_heat_injury():
#' hi <- predict_heat_injury(trace, wf, T_c = get_tcrit_summary(out2)$temp_median)
#' }
#' @export
extract_tdt <- function(workflow,
                        target_surv      = "relative",
                        t_ref            = 60,
                        TC_rate_range    = c(0.1, 1),
                        temp_grid        = NULL,
                        duration_grid    = NULL,
                        ndraws           = 1000,
                        time_multiplier  = NULL,
                        output_time_unit = "min",
                        lethal           = FALSE,
                        z_local          = FALSE,
                        seed             = NULL,
                        by               = NULL) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  if (length(TC_rate_range) != 2L ||
      any(!is.finite(TC_rate_range)) ||
      any(TC_rate_range <= 0) ||
      TC_rate_range[1] >= TC_rate_range[2])
    stop("TC_rate_range must be c(low, high) with 0 < low < high (% HI/hour).",
         call. = FALSE)

  # Resolve the model->output time multiplier from the workflow's duration_unit
  # when not supplied, so t_ref (in output units) is converted to model units
  # correctly regardless of whether the model was fit in minutes or hours.
  time_multiplier <- tdt_resolve_time_multiplier(time_multiplier, workflow$meta,
                                                 output_time_unit)
  ts   <- resolve_target_surv(target_surv)
  data <- workflow$data
  if (is.null(temp_grid)) {
    trange   <- range(data$temp, na.rm = TRUE)
    temp_grid <- seq(trange[1] - 2, trange[2] + 2, by = 0.05)
  }

  # Reproducibility: one set.seed() governs both the draw subsample (draw_ids,
  # matching the retired tdt_extract_pars selection) and the T_crit rate draws.
  if (!is.null(seed)) set.seed(seed)
  fit  <- get_brmsfit(workflow)
  Tbar <- workflow$meta$temp_mean
  bnd  <- workflow$meta$bounds
  by   <- tdt_resolve_by(workflow, by)
  did  <- tls_draw_ids(fit, ndraws)
  exposure_in_model_units <- t_ref / time_multiplier
  target <- log10(exposure_in_model_units)
  # The rate-multiplier T_crit is anchored at the 1-HOUR reference regardless of
  # t_ref (the -2.5z offset is derived from CTmax_1hr), so it needs its own target
  # at 60 min in model units. Matches tls() (R/tls.R: log10(60 / time_multiplier))
  # and equals `target` when t_ref = 60, so t_ref = 60 fits are unchanged.
  target_1hr <- log10(60 / time_multiplier)
  z_temp_grid <- sort(unique(data$temp)); z_temp_grid <- z_temp_grid[is.finite(z_temp_grid)]
  Lz <- length(z_temp_grid); Lc <- length(temp_grid); h <- 1e-3
  p  <- ts$prob %||% 0.5
  # Midpoint relative keeps the exact closed-form CTmax; direct/absolute use the
  # numerical inversion (the bent/backbone curve). z is always the central-diff
  # local z (exact for a linear midpoint), evaluated via posterior_linpred — so
  # there is no coefficient-name parsing and no parameterisation branch.
  mid_rel_closed <- ts$mode == "relative" &&
    !identical(workflow$meta$parameterization %||% "midpoint", "direct")

  if (isTRUE(lethal))
    message("extract_tdt(): T_crit reported under the rate-multiplier ",
            "definition; valid for damage-accumulation (lethal) endpoints ",
            "only. If your data are sublethal (knockdown, performance, ",
            "PSII, ...) ignore T_crit and supply T_c manually downstream.")

  # Single combined grid per group: temp_c = 0 (mid intercept) | z_temp_grid -/+ h
  # (local z) | temp_grid (CTmax inversion). One eval -> z and CTmax share draws,
  # so the T_crit pairing is genuinely joint.
  ztc      <- z_temp_grid - Tbar
  combo_tc <- c(0, ztc - h, ztc + h, temp_grid - Tbar)
  nd <- tls_build_grid(fit$data, by = by, temp = "temp_c", temp_grid = combo_tc)
  sp <- tls_eval_subpars(fit, nd, bnd, draw_ids = did, mode = ts$mode, p = p)

  # LT_x curve: descriptive intermediate (single-condition only; grouped fits use
  # the per-group z/CTmax/T_crit below, or tls()/predict_*(by=)). Computed HERE,
  # before the per-group T_crit rate draws, so derive_tdt_curve's posterior_linpred
  # subsample consumes RNG in the same order as the pre-refactor code and T_crit
  # stays reproducible against the existing fixtures.
  lt50_curve <- if (is.null(by)) derive_tdt_curve(
    workflow = workflow, temp_grid = temp_grid, duration_grid = duration_grid,
    target_surv = target_surv, ndraws = ndraws, time_multiplier = time_multiplier,
    output_time_unit = output_time_unit) else NULL

  per_group <- lapply(unique(nd$.grp), function(g) {
    gi      <- which(nd$.grp == g)              # 1 + 2Lz + Lc columns, in combo order
    c_mid0  <- gi[1]
    c_minus <- gi[1 + seq_len(Lz)]
    c_plus  <- gi[1 + Lz + seq_len(Lz)]
    c_ctmax <- gi[1 + 2L * Lz + seq_len(Lc)]

    z_obj <- tls_local_z(sp$logLT[, c_plus, drop = FALSE],
                         sp$logLT[, c_minus, drop = FALSE], h, z_temp_grid,
                         local = isTRUE(z_local))

    if (mid_rel_closed) {
      slope   <- rowMeans((sp$logLT[, c_plus, drop = FALSE] -
                           sp$logLT[, c_minus, drop = FALSE]) / (2 * h))
      Tc <- Tbar + (target - sp$logLT[, c_mid0]) / slope   # exact closed-form
    } else {
      Tc <- tls_invert_logLT(sp$logLT[, c_ctmax, drop = FALSE], target, temp_grid)
    }
    cd <- tibble::tibble(.draw = seq_along(Tc), temp = Tc) |> dplyr::filter(is.finite(temp))
    cq <- stats::quantile(cd$temp, c(0.025, 0.5, 0.975), na.rm = TRUE, names = FALSE)
    ctmax <- list(draws = cd,
                  summary = tibble::tibble(temp_lower = cq[1], temp_median = cq[2],
                                           temp_upper = cq[3]))

    tcb <- NULL
    if (isTRUE(lethal)) {
      # Anchor T_crit at the 1-hour CTmax (CTmax_1hr + z * log10(rate)), NOT the
      # CTmax at t_ref: the rate-multiplier offset is defined against CTmax_1hr,
      # so T_crit must be invariant to the reporting reference. extract_tdt
      # previously used CTmax at t_ref -> wrong T_crit for t_ref != 60. ct1 is
      # computed the same way as the CTmax above, just at the 1 h target.
      ct1 <- if (mid_rel_closed)
        Tbar + (target_1hr - sp$logLT[, c_mid0]) / slope
      else
        tls_invert_logLT(sp$logLT[, c_ctmax, drop = FALSE], target_1hr, temp_grid)
      cd1 <- tibble::tibble(.draw = seq_along(ct1), CTmax_1hr = ct1) |>
        dplyr::filter(is.finite(CTmax_1hr))
      paired <- dplyr::inner_join(dplyr::select(z_obj$draws, .draw, z), cd1,
                                  by = ".draw")
      paired$log10_rate <- stats::runif(nrow(paired), log10(TC_rate_range[1] / 100),
                                        log10(TC_rate_range[2] / 100))
      paired$T_crit <- paired$CTmax_1hr + paired$z * paired$log10_rate
      tcd <- tibble::tibble(.draw = paired$.draw, temp = paired$T_crit,
                            log10_rate = paired$log10_rate)
      tq  <- stats::quantile(tcd$temp, c(0.025, 0.5, 0.975), na.rm = TRUE, names = FALSE)
      tcb <- list(draws = tcd,
                  summary = tibble::tibble(TC_rate_low = TC_rate_range[1],
                                           TC_rate_high = TC_rate_range[2],
                                           temp_lower = tq[1], temp_median = tq[2],
                                           temp_upper = tq[3]))
    }
    list(gc = if (!is.null(by)) nd[gi[1], by, drop = FALSE] else NULL,
         z = z_obj, ctmax = ctmax, tcrit = tcb)
  })

  # Assemble. Single-condition (one group, by = NULL) -> the original nested
  # contract, byte-identical. Grouped -> each summary/draws gains moderator col(s).
  tag <- function(tb, gc) if (is.null(gc) || is.null(tb)) tb else cbind(gc, tb, row.names = NULL)
  z_summary       <- dplyr::bind_rows(lapply(per_group, function(pg) tag(pg$z$summary, pg$gc)))
  z_draws         <- dplyr::bind_rows(lapply(per_group, function(pg) tag(pg$z$draws, pg$gc)))
  ctmax_summary   <- dplyr::bind_rows(lapply(per_group, function(pg) tag(pg$ctmax$summary, pg$gc)))
  ctmax_draws     <- dplyr::bind_rows(lapply(per_group, function(pg) tag(pg$ctmax$draws, pg$gc)))
  z_local_block <- if (isTRUE(z_local)) list(
    draws   = dplyr::bind_rows(lapply(per_group, function(pg) tag(pg$z$local_draws, pg$gc))),
    summary = dplyr::bind_rows(lapply(per_group, function(pg) tag(pg$z$local_summary, pg$gc)))) else NULL
  t_crit_block <- if (isTRUE(lethal)) list(
    draws   = dplyr::bind_rows(lapply(per_group, function(pg) tag(pg$tcrit$draws, pg$gc))),
    summary = dplyr::bind_rows(lapply(per_group, function(pg) tag(pg$tcrit$summary, pg$gc)))) else NULL

  list(
    z          = list(draws = z_draws, summary = z_summary, local = z_local_block),
    CTmax      = list(draws = ctmax_draws, summary = ctmax_summary),
    T_crit     = t_crit_block,
    lt50_curve = lt50_curve,
    meta       = list(target_surv      = ts$label,
                      target_mode      = ts$mode,
                      target_prob      = ts$prob,
                      t_ref            = t_ref,
                      TC_rate_range    = TC_rate_range,
                      lethal           = isTRUE(lethal),
                      output_time_unit = output_time_unit,
                      by               = by)
  )
}
