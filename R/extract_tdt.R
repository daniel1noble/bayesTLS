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

# --- Shared posterior machinery for z and CTmax ------------------------------
# Both z and CTmax are deterministic transforms of the SAME population-level 4PL
# coefficient draws. extract_tdt() pulls ONE subsample with tdt_extract_pars()
# and feeds it to tdt_z_from_pars() and tdt_ctmax_from_pars(), so the two
# quantities share draws by construction — their correlation, and the joint
# pairing used for T_crit, are preserved with no draw-id bookkeeping and no
# random re-subsampling between quantities.

# Population-level 4PL coefficient draws, optionally subsampled to `ndraws`.
tdt_extract_pars <- function(workflow, ndraws = NULL) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  d <- posterior::as_draws_df(workflow$fit) |> as.data.frame()
  if (!is.null(ndraws) && is.finite(ndraws) && ndraws < nrow(d))
    d <- d[sort(sample.int(nrow(d), ndraws)), , drop = FALSE]
  d
}

# Build the `direct` reconstruction spec for a single-condition direct fit, or
# NULL for a midpoint fit. Grouped direct fits are redirected upstream by
# tdt_stop_if_grouped(); the `%in%` check is a belt-and-suspenders fallback so a
# caller that forgets the grouped guard still errors rather than silently reading
# the reference level. `p_fit` is the absolute backbone's anchoring probability
# baked in by make_4pl_formula (fit_4pl always anchors at 0.5).
tdt_direct_spec <- function(workflow, par_names) {
  if (!identical(workflow$meta$parameterization, "direct")) return(NULL)
  if (!all(c("b_CTmaxdev_Intercept", "b_logz_Intercept") %in% par_names))
    stop("This function supports single-condition direct fits; for grouped ",
         "direct fits use tls(object, by = ...).", call. = FALSE)
  list(ctmaxdev      = "b_CTmaxdev_Intercept",
       logz          = "b_logz_Intercept",
       fit_threshold = workflow$meta$threshold  %||% "relative",
       log10_tref    = workflow$meta$log10_tref %||% 0,
       p_fit         = workflow$meta$p          %||% 0.5)
}

# log10 LT_threshold(T) per draw, a [nrow(pars) x length(T)] matrix, from the
# coefficient draws. Relative threshold: log10 LT = mid(T). Absolute p:
# mid(T) + log((u - p)/(p - low))/k. A missing temp_c column (a dropped
# temp_effects slope) counts as zero slope — population level, no random effects.
tdt_loglt <- function(pars, Tbar, bnd, ts, T, direct = NULL) {
  np  <- nrow(pars)
  col <- function(nm) if (nm %in% names(pars)) pars[[nm]] else rep(0, np)
  lp  <- function(nlpar, Tv)
    outer(col(paste0("b_", nlpar, "_temp_c")), Tv - Tbar) +
      col(paste0("b_", nlpar, "_Intercept"))

  fit_abs <- !is.null(direct) && identical(direct$fit_threshold, "absolute")
  req_abs <- ts$mode != "relative"

  # Asymmetry correction C(T; p) = log((u - p)/(p - l))/k. TWO distinct
  # probabilities are involved and must NOT be conflated (a single Cmat at
  # p = ts$prob produced all-NA for an absolute-fit relative request, where
  # ts$prob is NA, and the wrong curve for a custom absolute request, where the
  # subtraction and addition cancelled):
  #   - p0 = the absolute backbone's anchoring probability baked in by
  #     make_4pl_formula (fit_4pl always uses 0.5). Used ONLY to UNDO the fit's
  #     baked-in midpoint correction (mid_curve = backbone - C(T; p0)).
  #   - p_req = ts$prob, the REQUESTED absolute level. Used ONLY to ADD the
  #     correction for an absolute request. They coincide only at p_req = 0.5.
  Cof <- NULL
  if (req_abs || fit_abs) {
    l <- bnd$low_min + bnd$low_w * stats::plogis(lp("lowraw", T))
    u <- bnd$up_min  + bnd$up_w  * stats::plogis(lp("upraw",  T))
    k <- exp(lp("logk", T))
    Cof <- function(pp) {
      arg <- (u - pp) / (pp - l)
      out <- log(arg) / k
      out[!is.finite(arg) | arg <= 0] <- NA_real_
      out
    }
  }

  # Curve midpoint mid(T). Midpoint parameterisation: the `mid` nlpar. Direct
  # parameterisation: reconstruct the linear backbone from CTmaxdev/logz
  # (= log10_tref - ((T - Tbar) - CTmaxdev)/exp(logz)); for an absolute FIT the
  # estimated midpoint is that backbone minus the fit-anchor correction C(T; p0).
  if (is.null(direct)) {
    mid <- lp("mid", T)
  } else {
    cd  <- col(direct$ctmaxdev)
    lz  <- col(direct$logz)
    mid <- direct$log10_tref - (outer(rep(1, np), T - Tbar) - cd) / exp(lz)
    if (fit_abs) mid <- mid - Cof(direct$p_fit %||% 0.5)
  }

  if (!req_abs) return(mid)
  mid + Cof(ts$prob)
}

# z per draw from coefficient draws. See derive_z() for the documented maths.
# `local = FALSE` skips the per-temperature z(T) breakdown (local_draws /
# local_summary); only the pooled z is needed unless the caller requests it.
tdt_z_from_pars <- function(pars, Tbar, bnd, ts, temp_grid,
                            probs = c(0.025, 0.5, 0.975), h = 1e-3,
                            local = TRUE, direct = NULL) {
  np <- nrow(pars)
  if (ts$mode == "relative" && is.null(direct)) {
    # Midpoint param: log10 LT_rel(T) = mid(T) is linear, z = -1 / b_mid_temp_c.
    zvec <- -1 / (if ("b_mid_temp_c" %in% names(pars)) pars$b_mid_temp_c
                  else rep(0, np))
    zloc <- matrix(zvec, nrow = np, ncol = length(temp_grid))
  } else {
    # Local slope of the (possibly bent / direct-reconstructed) log10 LT(T) by
    # central finite difference. Direct relative recovers z = exp(logz); direct
    # absolute / either-mode absolute reflects the curve's local slope.
    slope <- (tdt_loglt(pars, Tbar, bnd, ts, temp_grid + h, direct) -
              tdt_loglt(pars, Tbar, bnd, ts, temp_grid - h, direct)) / (2 * h)
    zloc  <- -1 / slope
    zloc[!is.finite(zloc)] <- NA_real_
  }
  z_pooled <- rowMeans(zloc, na.rm = TRUE)
  z_pooled[is.nan(z_pooled)] <- NA_real_

  draws <- tibble::tibble(.draw = seq_len(np), z = z_pooled) |>
    dplyr::filter(is.finite(z))
  qz <- stats::quantile(draws$z, probs, names = FALSE, na.rm = TRUE)
  # The per-temperature local-z breakdown is only surfaced by
  # extract_tdt(z_local = TRUE) and derive_z(); the pooled z above (rowMeans of
  # zloc) is unaffected by it. Building the np x nT tibble + group-by is the
  # dominant cost of this function, so skip it when it is not requested.
  if (isTRUE(local)) {
    local_draws <- tibble::tibble(
      .draw = rep(seq_len(np), times = length(temp_grid)),
      temp  = rep(temp_grid, each = np),
      z     = as.vector(zloc)
    ) |>
      dplyr::filter(is.finite(z))
    local_summary <- local_draws |>
      dplyr::group_by(temp) |>
      dplyr::summarise(
        z_median = stats::median(z),
        z_lower  = stats::quantile(z, probs[1], names = FALSE),
        z_upper  = stats::quantile(z, probs[3], names = FALSE),
        .groups  = "drop"
      )
  } else {
    local_draws   <- NULL
    local_summary <- NULL
  }
  list(draws = draws,
       summary = tibble::tibble(z_median = qz[2], z_lower = qz[1],
                                z_upper = qz[3]),
       local_draws = local_draws, local_summary = local_summary,
       target_surv = ts$label, temp_grid = temp_grid)
}

# CTmax per draw: temperature where log10 LT(T) = log10(exposure_model). Relative
# is the closed-form inverse of the linear mid(T); absolute finds the crossing of
# the closed-form LT curve over temp_grid by a vectorised inverse interpolation
# across draws, with an exact per-row fallback for non-monotone/non-finite curves.
tdt_ctmax_from_pars <- function(pars, Tbar, bnd, ts, exposure_model, temp_grid,
                                probs = c(0.025, 0.5, 0.975), direct = NULL) {
  np     <- nrow(pars)
  target <- log10(exposure_model)
  if (ts$mode == "relative" && is.null(direct)) {
    # Midpoint param: closed-form inverse of the linear mid(T).
    col <- function(nm) if (nm %in% names(pars)) pars[[nm]] else rep(0, np)
    Tc  <- Tbar + (target - col("b_mid_Intercept")) / col("b_mid_temp_c")
  } else {
    # log10 LT_p(T) is monotone-decreasing in T for the disjoint-bounds 4PL, so
    # it crosses `target` exactly once. Find that crossing by a vectorised
    # inverse linear interpolation across all draws at once: `j` is the column of
    # the last grid point still above target, and we interpolate within the
    # [j, j+1] interval. This replaces a per-draw stats::approx() loop and is
    # identical to it (verified to machine precision) for monotone curves.
    M  <- tdt_loglt(pars, Tbar, bnd, ts, temp_grid, direct)   # [np x nT]
    nT <- length(temp_grid)
    Tc <- rep(NA_real_, np)
    j  <- rowSums(M > target, na.rm = TRUE)
    ok <- is.finite(j) & j >= 1L & j <= (nT - 1L)
    iv <- which(ok); jj <- j[iv]
    y1 <- M[cbind(iv, jj)]; y2 <- M[cbind(iv, jj + 1L)]
    Tc[iv] <- temp_grid[jj] +
      (target - y1) * (temp_grid[jj + 1L] - temp_grid[jj]) / (y2 - y1)
    # Exact fallback: any draw whose curve is non-monotone or carries non-finite
    # values breaks the single-crossing assumption above, so reproduce the
    # original per-row interpolation for just those (rare) rows.
    dM  <- M[, -1, drop = FALSE] - M[, -nT, drop = FALSE]
    bad <- !is.finite(rowSums(M)) | rowSums(dM >= 0, na.rm = TRUE) > 0
    for (i in which(bad)) {
      y <- M[i, ]; fin <- is.finite(y)
      if (sum(fin) < 2L) { Tc[i] <- NA_real_; next }
      o <- order(y[fin])
      Tc[i] <- suppressWarnings(stats::approx(y[fin][o], temp_grid[fin][o],
                                              xout = target)$y)
    }
  }
  draws <- tibble::tibble(.draw = seq_len(np), temp = Tc) |>
    dplyr::filter(is.finite(temp))
  q <- stats::quantile(draws$temp, probs, names = FALSE, na.rm = TRUE)
  list(draws   = draws,
       summary = tibble::tibble(temp_lower = q[1], temp_median = q[2],
                                temp_upper = q[3]))
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
#' @return A list with `draws` (per-draw threshold temperatures; `target_surv`
#'         column is a character label), `summary` (quantile summary),
#'         `exposure_duration`, `target_surv` (the label).
#' @export
derive_temperature_for_duration <- function(workflow,
                                            exposure_duration,
                                            temp_grid,
                                            target_surv = "absolute",
                                            ndraws      = 1000,
                                            probs       = c(0.025, 0.5, 0.975),
                                            seed        = NULL) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  tdt_stop_if_grouped(workflow, "derive_temperature_for_duration()")

  ts     <- resolve_target_surv(target_surv)
  if (!is.null(seed)) set.seed(seed)   # reproducible draw subsample
  pars   <- tdt_extract_pars(workflow, ndraws)
  direct <- tdt_direct_spec(workflow, names(pars))
  ct     <- tdt_ctmax_from_pars(pars, workflow$meta$temp_mean,
                                workflow$meta$bounds, ts, exposure_duration,
                                temp_grid, probs, direct = direct)

  draws <- ct$draws |>
    dplyr::mutate(target_surv = ts$label) |>
    dplyr::select(.draw, target_surv, temp)
  summary <- ct$summary |>
    dplyr::mutate(target_surv = ts$label) |>
    dplyr::select(target_surv, temp_lower, temp_median, temp_upper)

  list(draws             = draws,
       summary           = summary,
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
#' @param h           Temperature step (°C) for the central finite difference
#'                    used in absolute mode. Default `1e-3`.
#' @param seed        Optional integer seeding the draw subsample (relevant only
#'                    when `ndraws` is set) for reproducibility. `NULL` (default)
#'                    leaves the RNG untouched.
#' @return A list with:
#'   - `draws`: tibble `(.draw, z)` — pooled per-draw z.
#'   - `summary`: tibble `(z_median, z_lower, z_upper)`.
#'   - `local_draws`: tibble `(.draw, temp, z)` — local z(T) per draw.
#'   - `local_summary`: tibble `(temp, z_median, z_lower, z_upper)`.
#'   - `target_surv`, `temp_grid`.
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
                     seed        = NULL) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  tdt_stop_if_grouped(workflow, "derive_z()")
  ts <- resolve_target_surv(target_surv)

  if (is.null(temp_grid)) temp_grid <- sort(unique(workflow$data$temp))
  temp_grid <- temp_grid[is.finite(temp_grid)]
  if (length(temp_grid) < 1L) stop("temp_grid is empty.", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)   # reproducible draw subsample (when ndraws set)
  pars   <- tdt_extract_pars(workflow, ndraws)
  direct <- tdt_direct_spec(workflow, names(pars))
  tdt_z_from_pars(pars, workflow$meta$temp_mean, workflow$meta$bounds,
                  ts, temp_grid, probs, h, direct = direct)
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
                        seed             = NULL) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  # Grouped fits (CTmax/z by a moderator) -> per-group derivation via tls();
  # detected coding-independently so `~ 0 + G` and `~ 1 + G` both redirect.
  tdt_stop_if_grouped(workflow, "extract_tdt()")
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

  # Validate up front so both helpers receive a normalised label.
  ts <- resolve_target_surv(target_surv)

  data <- workflow$data
  if (is.null(temp_grid)) {
    trange   <- range(data$temp, na.rm = TRUE)
    temp_grid <- seq(trange[1] - 2, trange[2] + 2, by = 0.05)
  }

  # Extract the population-level 4PL coefficient posterior ONCE and subsample
  # ONCE, then compute BOTH z and CTmax from this single set of draws. They
  # therefore share draws by construction — their correlation, and the per-
  # `.draw` pairing used for T_crit below, are preserved without any draw-id
  # bookkeeping. z is read directly (no regression): relative => -1/b_mid_temp_c
  # per draw; absolute => pooled local z over the observed assay temperatures.
  # Reproducibility: seed the draw subsample (below) and the T_crit rate draws
  # (further down) from one stream, so a supplied seed fixes the whole result.
  if (!is.null(seed)) set.seed(seed)
  pars        <- tdt_extract_pars(workflow, ndraws)
  Tbar        <- workflow$meta$temp_mean
  bnd         <- workflow$meta$bounds

  # Direct CTmax/z parameterisation: there is no `mid` coefficient, so mid(T) is
  # reconstructed from the CTmaxdev/logz coefficient draws (grouped fits already
  # redirected by tdt_stop_if_grouped above; midpoint fits get NULL).
  direct <- tdt_direct_spec(workflow, names(pars))

  z_temp_grid <- sort(unique(data$temp))
  z_temp_grid <- z_temp_grid[is.finite(z_temp_grid)]

  z_obj <- tdt_z_from_pars(pars, Tbar, bnd, ts, z_temp_grid,
                           local = isTRUE(z_local), direct = direct)

  # Express the t_ref reference duration back in model time units so the
  # inverse-4PL lookup uses the same scale as the fitted model.
  exposure_in_model_units <- t_ref / time_multiplier
  ctmax <- tdt_ctmax_from_pars(pars, Tbar, bnd, ts, exposure_in_model_units,
                               temp_grid, direct = direct)

  # LT_x curve retained as descriptive output (plotting / inspection); it is no
  # longer used to compute z or CTmax. The threshold is set by target_surv;
  # `"relative"` (default) returns mid(T) directly.
  lt50_curve <- derive_tdt_curve(
    workflow         = workflow,
    temp_grid        = temp_grid,
    duration_grid    = duration_grid,
    target_surv      = target_surv,
    ndraws           = ndraws,
    time_multiplier  = time_multiplier,
    output_time_unit = output_time_unit
  )

  t_crit_block <- NULL
  if (isTRUE(lethal)) {
    message("extract_tdt(): T_crit reported under the rate-multiplier ",
            "definition; valid for damage-accumulation (lethal) endpoints ",
            "only. If your data are sublethal (knockdown, performance, ",
            "PSII, ...) ignore T_crit and supply T_c manually downstream.")

    # T_crit via rate-multiplier integration. For each posterior draw, sample
    # r* uniformly on log10 across TC_rate_range, then compute
    # T_crit = CTmax + z * log10(r*/100). z and CTmax were computed from the
    # same `pars` draws, so pairing by `.draw` is genuinely joint.
    z_df     <- z_obj$draws  |> dplyr::select(.draw, z)
    ctmax_df <- ctmax$draws  |> dplyr::select(.draw, CTmax_temp = temp)
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
    z          = list(draws   = z_obj$draws,
                      summary = z_obj$summary,
                      local   = if (isTRUE(z_local)) {
                        list(draws = z_obj$local_draws,
                             summary = z_obj$local_summary)
                      } else {
                        NULL
                      }),
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
