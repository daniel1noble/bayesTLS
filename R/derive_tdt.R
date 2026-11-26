# Classical TDT quantities (z, CTmax, T_crit) derived from a fitted 4PL. The LT
# at a survival threshold is read per draw from the shared engine
# (tls_eval_subpars()): the 4PL `mid` for the relative threshold (default), and
# `mid` plus the model's closed-form offset for an absolute survival
# probability. Temperatures at a fixed exposure use the exact closed form when
# that LT curve is a straight line in temperature, and invert it on a grid
# otherwise.

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
#'   numerical inversion. When `low` is about 0 and `up` about 1 it coincides with
#'   the classical LT50.
#' - `"absolute"` (or a numeric `p` in `(0, 1)`): the duration at which
#'   survival crosses the literal probability `p` (0.5 by default), read
#'   exactly from each draw's fitted curve (the midpoint plus the model's own
#'   closed-form correction, as in [tls()]), so no duration grid bounds the
#'   answer. A draw whose curve never reaches `p` at that temperature (`p`
#'   outside its asymptotes) has no crossing; such draws are dropped from the
#'   summary with a warning.
#'
#' This is the **horizontal** read of the survival surface: fix a survival
#' threshold, read off the time required to reach it at each temperature.
#'
#' @param workflow         Fitted `bayes_tls`.
#' @param temp_grid        Numeric vector of temperatures (°C).
#' @param duration_grid    Deprecated and ignored, with a message if supplied
#'                         (in either threshold mode): absolute crossings are
#'                         now computed exactly per draw, so no search grid is
#'                         needed. It used to bound the numerical search,
#'                         silently dropping draws whose crossing lay beyond
#'                         it.
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
#' @param by               Optional moderator column(s) for per-group LT curves.
#'                         `NULL` (default) uses the fit's recorded moderators;
#'                         a single-condition fit returns one ungrouped curve.
#' @return A list with `draws` (per-draw threshold durations; `target_surv`
#'         column is a character label; a duration too long or too short to
#'         represent is `Inf` or 0), `summary` (quantile summary by
#'         temperature, plus moderator columns for grouped fits), `target_surv`
#'         (the label), `time_multiplier`, `output_time_unit`, and `by`.
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
                             output_time_unit = "min",
                             by               = NULL) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  tdt_check_ndraws(ndraws)
  temp_grid <- tdt_finite_grid(temp_grid)
  if (!is.null(duration_grid))
    message("`duration_grid` is no longer used: derive_tdt_curve() reads ",
            "absolute crossings exactly from each draw's fitted curve.")

  # brms::posterior_linpred() errors if `ndraws` exceeds the posterior size, so
  # clamp it for the relative branch's single posterior_linpred(ndraws =) call
  # (a fit with fewer draws, e.g. 2 chains x 400, would otherwise crash). The
  # absolute branch subsamples through tls_draw_ids(), which uses all draws
  # when `ndraws` is at least the posterior size.
  if (!is.null(ndraws)) {
    total <- tryCatch(brms::ndraws(workflow$fit), error = function(e) NA_integer_)
    if (is.finite(total)) ndraws <- min(ndraws, total)
  }

  time_multiplier <- tdt_resolve_time_multiplier(time_multiplier, workflow$meta,
                                                 output_time_unit)
  ts <- resolve_target_surv(target_surv)
  by <- tdt_resolve_by(workflow, by)

  if (ts$mode == "relative") {
    # Direct shortcut: log10(t_relative) = mid(T) per draw. No grid search.
    nd <- new_tdt_grid(workflow, temps = temp_grid, durations = 1, by = by)
    pp_mid <- brms::posterior_linpred(workflow$fit, newdata = nd,
                                       nlpar = "mid", re_formula = NA,
                                       ndraws = ndraws)
    # pp_mid is [ndraws x nrow(nd)] of log10(t) in model time units.
    duration_model_mat <- 10 ^ pp_mid

    draw_list <- vector("list", nrow(nd))
    for (i in seq_len(nrow(nd))) {
      t_i  <- nd$temp[i]
      dmod <- duration_model_mat[, i]
      d <- data.frame(
        .draw            = seq_along(dmod),
        temp             = t_i,
        target_surv      = ts$label,
        duration_model   = dmod,
        duration_out     = dmod * time_multiplier,
        stringsAsFactors = FALSE
      )
      if (!is.null(by)) {
        d <- cbind(nd[rep(i, nrow(d)), by, drop = FALSE], d, row.names = NULL)
      }
      draw_list[[i]] <- d
    }
  } else {
    # Exact per-draw LT from the shared engine (midpoint + the model's
    # closed-form offset), the same logLT that tls(), derive_z() and
    # derive_temperature_for_duration() read. The old numerical search on a
    # duration grid silently dropped every draw whose crossing lay beyond the
    # grid, so cool temperatures were summarised from a biased few draws.
    fit <- workflow$fit
    nd  <- new_tdt_grid(workflow, temps = temp_grid, durations = 1, by = by)
    did <- tls_draw_ids(fit, ndraws)   # one subsample for all four sub-parameters
    sp  <- tls_eval_subpars(fit, nd, tdt_workflow_bounds(workflow), draw_ids = did,
                            mode = "absolute", p = ts$prob)
    duration_model_mat <- 10 ^ sp$logLT   # NA where p lies outside the asymptotes

    draw_list <- vector("list", nrow(nd))
    for (i in seq_len(nrow(nd))) {
      dmod <- duration_model_mat[, i]
      d <- data.frame(
        .draw            = seq_along(dmod),
        temp             = nd$temp[i],
        target_surv      = ts$label,
        duration_model   = dmod,
        duration_out     = dmod * time_multiplier,
        stringsAsFactors = FALSE
      )
      if (!is.null(by)) {
        d <- cbind(nd[rep(i, nrow(d)), by, drop = FALSE], d, row.names = NULL)
      }
      draw_list[[i]] <- d
    }
  }

  # Keep durations too long or too short to represent (10^logLT overflows to
  # Inf, e.g. a target just above the lower asymptote, or underflows to 0):
  # dropping them would silently cut the most extreme draws from the summary.
  # Only undefined draws (NA: target outside the asymptotes) are dropped.
  draws <- dplyr::bind_rows(draw_list) |>
    dplyr::filter(!is.na(duration_model))

  group_vars <- c(if (!is.null(by)) by, "target_surv", "temp")
  summary <- draws |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
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
       output_time_unit = output_time_unit,
       by               = by)
}

#' Temperature at which survival equals a target after a fixed exposure
#'
#' The **vertical** read of the survival surface: fix an exposure duration,
#' find the temperature at which the posterior survival reaches the chosen
#' threshold. Returns one temperature per posterior draw.
#'
#' Threshold modes (via `target_surv`) match [derive_tdt_curve()]:
#'
#' - `"relative"` (default) → temperature at which
#'   `mid(T) = log10(exposure_duration)` per draw.
#' - `"absolute"` (= 0.5) or numeric `p` in `(0, 1)` → temperature at which
#'   the LT at the literal survival probability `p` (the midpoint plus the
#'   model's closed-form offset) equals the exposure.
#'
#' Each draw's temperature comes from the exact closed form whenever its
#' `log10(LT)` is a straight line in temperature: every relative threshold of
#' a fit whose midpoint is linear in temperature (midpoint and direct fits
#' alike), and an absolute threshold when the asymptotes and steepness do not
#' change with temperature. A curve that bends (temperature effects on the
#' asymptotes or steepness, under an absolute threshold) is inverted on
#' `temp_grid`, refined as in [tls()] (each interval split into equal steps of
#' at most 0.05 degrees C, but into no more than 20, so an interval wider than
#' 1 degree C is refined to a twentieth of its width), by linear interpolation
#' between the refined points, taking the crossing
#' where the LT shortens through the exposure as temperature rises; draws that
#' do not cross on it are searched outward beyond each end, and draws left
#' without a single crossing are reported in a warning. Where an absolute
#' target lies outside the fitted asymptotes at a temperature, the LT there is
#' treated as zero or infinite, so a crossing next to it is still found.
#' [tls()] reads `CTmax` at `t_ref` the same way.
#'
#' @param workflow         Fitted `bayes_tls`.
#' @param exposure_duration Numeric scalar — the fixed duration (model units).
#' @param temp_grid        Numeric vector of temperatures (°C) to search over
#'                         when the curve bends; it is extended beyond each end
#'                         by the data's temperature range (or its own width,
#'                         if wider; a single temperature is extended at a
#'                         fine step). Not needed for a straight curve, which
#'                         is solved exactly. `NULL` (default): 11 points over
#'                         the observed temperatures.
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
                                            temp_grid   = NULL,
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

  # Evaluate logLT per group at temp_c = 0 and 1 (the closed-form anchor and
  # slope), at a probe temp_c = 0.5 and over the search grid (the inversion of
  # a bent curve). The probe gives the straightness check a third temperature
  # even when the grid is only temp_c = 0 and/or 1, so such a grid cannot make
  # a bent curve look straight. Rows within a group keep this order
  # (tls_build_grid()), so columns 1-3 are temp_c = 0, 1 and 0.5 and the rest
  # the grid.
  if (is.null(temp_grid))
    temp_grid <- seq(min(workflow$data$temp), max(workflow$data$temp), length.out = 11)
  tg <- sort(unique(tdt_finite_grid(temp_grid) - tbar))
  bounds <- tdt_workflow_bounds(workflow)
  pp     <- if (is.na(ts$prob)) 0.5 else ts$prob
  width  <- diff(range(fit$data$temp_c, na.rm = TRUE))
  nd <- tls_build_grid(fit$data, by = by, temp = "temp_c", temp_grid = c(0, 1, 0.5, tg))
  sp <- tls_eval_subpars(fit, nd, bounds, draw_ids = did, mode = ts$mode, p = pp)

  # Per group: a straight logLT(T) is solved exactly at any temperature; a
  # bent one is inverted on the grid (sign-safely at an absolute threshold),
  # then the draws without a crossing there are searched outward, all groups
  # together (tls_extend_crossings()).
  p_abs  <- if (identical(ts$mode, "absolute")) pp else NA_real_
  model  <- tdt_model(fit)
  groups <- unique(nd$.grp)
  Tcs <- list(); bent <- list()
  for (g in groups) {
    gi <- which(nd$.grp == g)
    M  <- sp$logLT[, gi, drop = FALSE]
    if (tls_group_is_linear(M, c(0, 1, 0.5, tg))) {
      Tcs[[g]] <- tbar + (target - M[, 1]) / (M[, 2] - M[, 1])
    } else {
      spg <- lapply(sp[c("low", "up", "k", "mid", "logLT")],
                    function(m) m[, gi[-(1:3)], drop = FALSE])
      R   <- if (is.na(p_abs)) NULL else tls_surv_resid(spg, target, p_abs, model)
      Tcs[[g]] <- tls_invert_logLT(tls_clamp_logLT(spg, p_abs), target, tg + tbar, R = R)
      bent[[g]] <- list(newrow = nd[gi[1], , drop = FALSE], tg = tg, sp = spg,
                        crossings = list(v = Tcs[[g]]), targets = list(v = target))
    }
  }
  uncrossed <- character(0)
  if (length(bent)) {
    bent <- tls_extend_crossings(fit, bent, "temp_c", width = width, Tbar = tbar,
                                 draw_ids = did, bounds = bounds, mode = ts$mode,
                                 p = p_abs, model = model)
    for (g in names(bent)) {
      Tcs[[g]] <- bent[[g]]$crossings$v
      n_miss <- sum(!is.finite(Tcs[[g]]))
      if (n_miss > 0L)
        uncrossed <- c(uncrossed, sprintf("%s%d of %d draws",
                                          if (is.null(by)) "" else paste0(g, ": "),
                                          n_miss, length(Tcs[[g]])))
    }
  }
  per_group <- lapply(groups, function(g) {
    gi <- which(nd$.grp == g)
    Tc <- Tcs[[g]]
    d <- tibble::tibble(.draw = seq_along(Tc), target_surv = ts$label, temp = Tc) |>
      dplyr::filter(is.finite(temp))
    q <- stats::quantile(d$temp, probs, names = FALSE, na.rm = TRUE)
    s <- tibble::tibble(target_surv = ts$label, temp_lower = q[1],
                        temp_median = q[2], temp_upper = q[3])
    if (!is.null(by)) {
      gc <- nd[gi[1], by, drop = FALSE]
      d <- cbind(gc[rep(1L, nrow(d)), , drop = FALSE], d, row.names = NULL)  # 0 rows if none crossed
      s <- cbind(gc, s, row.names = NULL)
    }
    list(draws = d, summary = s)
  })
  if (length(uncrossed)) {
    ext <- range(tls_extend_grid(tg, width)) + tbar
    warning(sprintf(paste0(
      "derive_temperature_for_duration(): some draws have no single crossing ",
      "of the exposure within %.4g to %.4g degrees C (temp_grid extended beyond ",
      "each end by the data's temperature range or its own width, whichever is ",
      "larger): %s. They are left out of the summary, which can then ",
      "under-represent the extremes. A draw has none when its curve stays on ",
      "one side of the target over the search (the crossing lies further out, ",
      "or the target is outside the fitted asymptotes there) or crosses it more ",
      "than once; a wider temp_grid helps only in the first case."),
      ext[1], ext[2], paste(uncrossed, collapse = "; ")), call. = FALSE)
  }

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
#' - **Absolute threshold** \eqn{p}: the LT curve gains the fitted model's
#'   offset from the midpoint, \eqn{\log_{10}\mathrm{LT}_p(T) = \mathrm{mid}(T) +
#'   C_p(T)}. For the log-logistic ([fit_4pl()])
#'   \eqn{C_p = \tfrac{1}{k(T)}\log\tfrac{u(T)-p}{p-\ell(T)}}; the log-normal
#'   and Weibull models ([fit_lognormal()], [fit_weibull()]) use their own
#'   closed forms. When \eqn{\ell}, \eqn{u} or
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
  sp <- tls_eval_subpars(fit, nd, tdt_workflow_bounds(workflow), draw_ids = did,
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
