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
#' @param ndraws           Posterior draws to use; `NULL` (default) uses the full
#'                       posterior. Subsampling can shift results slightly off
#'                       the model, so the default keeps derived quantities
#'                       consistent with the fit; pass an integer for speed.
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
#' @examples
#' \dontrun{
#' wf  <- fit_4pl(std)
#' crv <- derive_tdt_curve(wf, temp_grid = c(38, 40, 42))  # relative LT curve
#' crv$summary
#' }
#' @export
derive_tdt_curve <- function(workflow,
                             temp_grid,
                             duration_grid    = NULL,
                             target_surv      = "relative",
                             ndraws           = NULL,
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
#' This is the primitive used by [tls()] to derive CTmax at `t_ref`.
#'
#' @param workflow         Fitted `bayes_tls`.
#' @param exposure_duration Numeric scalar — the fixed duration (model units).
#' @param temp_grid        Numeric vector of temperatures to search over.
#' @param target_surv      Threshold mode. `"relative"` (default), `"absolute"`,
#'                         or a numeric in `(0, 1)`.
#' @param ndraws           Posterior draws to use; `NULL` (default) uses the
#'                         full posterior. Pass an integer to subsample.
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
#' @examples
#' \dontrun{
#' wf <- fit_4pl(std)
#' # Temperature giving the relative midpoint threshold after a 60-unit exposure:
#' tt <- derive_temperature_for_duration(wf, exposure_duration = 60,
#'                                        temp_grid = seq(36, 44, by = 0.1))
#' tt$summary
#' }
#' @export
derive_temperature_for_duration <- function(workflow,
                                            exposure_duration,
                                            temp_grid,
                                            target_surv = "relative",
                                            ndraws      = NULL,
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
