# Generic, moderator-aware TDT extraction for arbitrary fitted 4PL models.
# Unlike extract_tdt() (which targets the standard fit_4pl() workflow), tls()
# evaluates the four 4PL sub-parameters at a user-defined moderator x
# temperature grid via brms::posterior_linpred(nlpar = ), so it works on
# hand-written brms models with moderators on ANY sub-parameter and arbitrary
# random-effect structures. No coefficient names are parsed.

#' Thermal load sensitivity summaries from any fitted 4PL model
#'
#' One call that derives the classical TDT quantities — thermal sensitivity
#' `z`, `CTmax`, and (for lethal endpoints) `T_crit` — **per moderator group**
#' from a fitted joint 4PL, including hand-written `brms` models with moderators
#' (sex, life stage, clone, ...) on any sub-parameter. The four sub-parameters
#' are evaluated at a moderator x temperature grid with
#' [brms::posterior_linpred()], then `z`, `CTmax` and `T_crit` are derived from
#' the same posterior draws (so they are mutually consistent and `T_crit` pairs
#' `CTmax` and `z` per draw).
#'
#' @param object A fitted `bayes_tls` workflow or a `brmsfit` 4PL whose
#'   non-linear parameters are the (reparameterised) asymptotes, steepness and
#'   midpoint.
#' @param by Character vector of moderator columns defining the groups reported
#'   separately (e.g. `"sex"`). `NULL` (default) pools to a single group.
#' @param params Quantities to return: `"all"` (z, CTmax, and T_crit when
#'   `lethal = TRUE`), or a subset of `c("z", "ctmax", "tcrit")`.
#' @param target_surv Survival threshold the LT curve is read at, shared with
#'   [extract_tdt()]: `"relative"` (default; the per-draw midpoint
#'   `(low + up)/2`), `"absolute"` (literal LT50), or a numeric `p` in `(0, 1)`
#'   for an absolute LT`p`. The default is `"relative"` so `tls()` and
#'   [extract_tdt()] return the same z/CTmax for a given fit.
#' @param mode,p **Deprecated**, superseded by `target_surv`. If supplied,
#'   `mode = "relative"` maps to `target_surv = "relative"` and
#'   `mode = "absolute"` (with `p`) to `target_surv = p`. Kept for back-compat.
#' @param t_ref Reference exposure duration for `CTmax`, in minutes. If `NULL`
#'   (default), inherit the fit's own reference exposure (`meta$t_ref`) and
#'   report it via a message; a raw `brmsfit` with no recorded reference
#'   falls back to 60.
#'   Converted to the model's time scale via `time_multiplier`.
#' @param time_multiplier Multiplier from the model's time unit to minutes.
#'   `NULL` (default) derives it from a `bayes_tls` workflow's `duration_unit`
#'   (e.g. 60 for an hours model); a raw `brmsfit` with no metadata defaults to
#'   1, so `t_ref` is then in the model's own time units.
#' @param lethal Logical; with `TRUE`, T_crit (rate-multiplier) is available.
#'   T_crit is meaningful only for damage-accumulation (lethal) endpoints.
#' @param TC_rate_range Length-2 HI-rate floor range (% per hour) for T_crit.
#' @param temp Name of the centred temperature column. Default `"temp_c"`.
#' @param temp_mean Centring constant mapping `temp` back to temperature, needed
#'   for `CTmax`/`T_crit`. Taken from a `bayes_tls` workflow's metadata when
#'   available; supply it for a raw `brmsfit`.
#' @param temp_grid Temperatures (on the `temp` scale) at which the LT curve is
#'   evaluated: `z` is the local `z(T)` pooled over this grid and `CTmax` is the
#'   per-draw crossing of `t_ref` inverted on it. Default: 11 points over the
#'   observed range. For an absolute threshold, ensure the grid brackets the
#'   `CTmax` crossing (extend it if `CTmax` falls outside the observed range).
#' @param re_formula Passed to [brms::posterior_linpred()]; default `NA`
#'   (population level). Include group-level terms (e.g. clone) for per-group draws.
#' @param lower,upper Response bounds of the disjoint-bounds reparameterisation.
#' @param nlpars Non-linear parameter names, ordered low / up / k / mid (raw scale).
#' @param ndraws Optional number of posterior draws to subsample.
#' @param newdata Optional explicit moderator x temperature grid; overrides the
#'   `by`/`temp_grid` construction.
#' @param probs Summary quantiles (lower, median, upper). Default `c(.025,.5,.975)`.
#' @param seed Optional integer seeding the draw subsample (`ndraws`) and the
#'   T_crit rate draws for reproducibility. `NULL` (default) leaves the RNG alone.
#' @param ... Additional arguments passed on to [tls()] (used by the
#'   `tls_z()`, `tls_ctmax()` and `tls_tcrit()` convenience wrappers).
#' @return A `tls` object: `$summary` (per-group, per-quantity median + interval),
#'   `$draws` (per-group, per-quantity posterior draws), and `$meta`.
#' @details
#' **Shared engine with [extract_tdt()].** `tls()` and [extract_tdt()] derive
#' `z`, `CTmax` and `T_crit` through the same internal engine, so on a common
#' temperature grid they agree to numerical tolerance for every threshold. `z`
#' is the central-difference local `z(T)` (negative reciprocal of the local
#' slope of `log10(LT)` on temperature) pooled over `temp_grid`; `CTmax` is the
#' per-draw temperature at which the fitted LT curve crosses `t_ref`, from the
#' exact closed form when the relative midpoint is linear and from a true
#' per-draw inversion of the (possibly bent) curve otherwise. Under
#' `target_surv = "relative"` (the default) `log10(LT) = mid(T)` is linear, so
#' `z` is constant in temperature and equals `-1 / b_mid_temp_c`; under an
#' absolute threshold with temperature effects on the asymptotes or steepness
#' the curve bends and both the local `z(T)` and the inversion account for it,
#' with no linear approximation. The two functions differ only in defaults, not
#' engine: `tls()` summarises `z` and `CTmax` on the single `temp_grid`, whereas
#' [extract_tdt()] pools `z` over the observed assay temperatures and inverts
#' `CTmax` on a finer extended grid.
#' @examples
#' \dontrun{
#' tls(joint_sex_fit, by = "sex", lethal = TRUE, temp_mean = 36.1)  # z, CTmax, T_crit per sex
#' tls(wf_leaf, params = "z")                                       # z only, workflow
#' }
#' @export
tls <- function(object, by = NULL, params = "all",
                target_surv = "relative",
                t_ref = NULL, time_multiplier = NULL,
                lethal = FALSE, TC_rate_range = c(0.1, 1),
                temp = "temp_c", temp_mean = NULL, temp_grid = NULL,
                re_formula = NA, lower = 0, upper = 1,
                nlpars = c("lowraw", "upraw", "logk", "mid"),
                ndraws = NULL, newdata = NULL,
                probs = c(0.025, 0.5, 0.975), seed = NULL,
                mode = NULL, p = NULL) {
  # `mode`/`p` are superseded by `target_surv` (the same argument extract_tdt()
  # uses); map them when supplied so existing calls keep working. target_surv
  # defaults to "relative" to MATCH extract_tdt() -- tls() previously defaulted
  # to absolute, so the two engines returned different headline z/CTmax for the
  # same fit (the documented inconsistency this unifies).
  if (!is.null(mode)) {
    mode <- match.arg(mode, c("absolute", "relative"))
    target_surv <- if (identical(mode, "relative")) "relative" else (p %||% 0.5)
  } else if (!is.null(p)) {
    target_surv <- p
  }
  ts   <- resolve_target_surv(target_surv)
  mode <- ts$mode
  p    <- if (is.na(ts$prob)) 0.5 else ts$prob

  # Cheap argument validation first (no fit needed).
  if (identical(params, "all")) params <- c("z", "ctmax", if (lethal) "tcrit")
  params <- match.arg(params, c("z", "ctmax", "tcrit"), several.ok = TRUE)
  if ("tcrit" %in% params && !lethal) {
    stop("T_crit is only defined for lethal endpoints; set `lethal = TRUE`. ",
         "The rate-multiplier T_crit is not meaningful for performance endpoints.",
         call. = FALSE)
  }

  if (inherits(object, "bayes_tls")) {
    fit  <- get_brmsfit(object)
    meta <- object$meta
    if (is.null(temp_mean)) temp_mean <- meta$temp_mean
  } else if (inherits(object, "brmsfit")) {
    fit  <- object
    meta <- list()
  } else {
    stop("`object` must be a bayes_tls workflow or a brmsfit.", call. = FALSE)
  }
  # Inherit the fit's own reference exposure when t_ref is omitted (and announce
  # it), so a fit done at e.g. 4 h is not silently read at the 60-min default.
  t_ref <- tdt_resolve_t_ref(t_ref, meta)
  # t_ref is in output minutes; convert to the model's time scale (raw fits with
  # no duration_unit metadata default to time_multiplier = 1, i.e. model units).
  time_multiplier <- tdt_resolve_time_multiplier(time_multiplier, meta, "min")

  if (any(c("ctmax", "tcrit") %in% params) && is.null(temp_mean)) {
    stop("`temp_mean` is required for CTmax / T_crit (the centring constant that ",
         "maps `", temp, "` back to temperature). For a bayes_tls workflow it is ",
         "read from metadata.", call. = FALSE)
  }

  mdata <- fit$data
  if (!temp %in% names(mdata)) {
    stop("Temperature column `", temp, "` not found in the model data.",
         call. = FALSE)
  }

  # --- build the moderator x temperature grid (shared engine) ----------------
  base_grid <- tls_build_grid(mdata, by = by, temp = temp,
                              temp_grid = temp_grid, newdata = newdata)

  # tls() derives z / CTmax / T_crit through the SAME engine as extract_tdt()
  # (tls_zct_draws): the central-difference local z pooled over the grid, and
  # either the exact closed form (a linear relative midpoint) or a true per-draw
  # inversion of the -- possibly bent -- absolute LT curve. To do that we
  # evaluate logLT at, per group, temp_c = 0 (the closed-form midpoint anchor),
  # the grid +/- h (local z), and the grid itself (the CTmax inversion). h is the
  # same finite-difference step extract_tdt() uses. This replaces the earlier
  # single global least-squares slope, which was exact only for a linear
  # (relative) LT curve and a ~linear approximation of a bent (absolute) one.
  h          <- 1e-3
  log_tref   <- log10(t_ref / time_multiplier)        # CTmax target (model units)
  target_1hr <- log10(60   / time_multiplier)         # T_crit 1-hour anchor
  # A linear relative midpoint has a closed-form crossing; an absolute threshold
  # or a `direct` fit is inverted numerically. Mirrors extract_tdt()'s switch.
  mid_rel_closed <- identical(mode, "relative") &&
    !identical(meta$parameterization %||% "midpoint", "direct")

  groups   <- unique(base_grid$.grp)
  grp_cols <- list(); aug_parts <- list(); offset <- 0L
  for (g in groups) {
    rows   <- which(base_grid$.grp == g)
    tg     <- sort(unique(base_grid[[temp]][rows]))   # centred grid temperatures
    L      <- length(tg)
    modrow <- base_grid[rows[1], , drop = FALSE]      # moderators + filled columns
    aug_tc <- c(0, tg - h, tg + h, tg)                # mid0 | minus | plus | ctmax
    a         <- modrow[rep(1L, length(aug_tc)), , drop = FALSE]
    a[[temp]] <- aug_tc
    a$.grp    <- g
    aug_parts[[length(aug_parts) + 1L]] <- a
    grp_cols[[g]] <- list(
      tg = tg, L = L,
      mod   = if (is.null(by)) NULL else base_grid[rows[1], by, drop = FALSE],
      mid0  = offset + 1L,
      minus = offset + 1L + seq_len(L),
      plus  = offset + 1L + L + seq_len(L),
      ctmax = offset + 1L + 2L * L + seq_len(L))
    offset <- offset + 1L + 3L * L
  }
  aug <- do.call(rbind, aug_parts)
  rownames(aug) <- NULL

  # Reproducibility: seed the posterior_linpred draw subsample (when `ndraws` is
  # set) and the T_crit rate draws below from one stream. The augmented grid adds
  # rows, not draws, so the subsample -- and hence seeded T_crit -- is unchanged.
  if (!is.null(seed)) set.seed(seed)

  # --- evaluate sub-parameters at every grid row (shared engine) -------------
  logLT <- tls_eval_subpars(fit, aug, compute_4pl_bounds(lower, upper),
                            nlpars = nlpars, re_formula = re_formula,
                            ndraws = ndraws, mode = mode, p = p)$logLT
  np <- nrow(logLT)

  want_z     <- "z"     %in% params
  want_ctmax <- "ctmax" %in% params
  want_tcrit <- "tcrit" %in% params

  # --- derive per-group quantities from the shared draws ---------------------
  summ <- list(); drw <- list()
  allna_groups <- character(0)   # groups whose z/CTmax draws are all non-finite
  for (g in groups) {
    gi    <- grp_cols[[g]]
    gcols <- gi$mod
    add <- function(q, v) {
      # Summarise on the finite draws only. A near-flat LT curve sends
      # z = -1/slope to +/-Inf, a curve that never crosses t_ref within the grid
      # gives an NA CTmax, and a single-temperature group is NA'd below; without
      # this, stats::quantile() aborts ("missing values not allowed") and takes
      # down the whole tls() call. The raw draws (Inf/NaN included) are still
      # returned in `$draws`.
      vf <- v[is.finite(v)]
      if (!length(vf)) allna_groups <<- union(allna_groups, as.character(g))
      qs <- if (length(vf))
              stats::quantile(vf, probs[c(1, 3)], names = FALSE) else c(NA_real_, NA_real_)
      s <- data.frame(quantity = q,
                      median = if (length(vf)) stats::median(vf) else NA_real_,
                      lower = qs[1], upper = qs[2], row.names = NULL)
      d <- data.frame(quantity = q, .draw = seq_along(v), value = v,
                      row.names = NULL)
      if (!is.null(gcols)) {                       # data.frame recycles the 1-row gcols
        s <- data.frame(gcols, s, row.names = NULL)
        d <- data.frame(gcols, d, row.names = NULL)
      }
      summ[[length(summ) + 1L]] <<- s
      drw[[length(drw) + 1L]]   <<- d
    }

    if (gi$L < 2L) {
      # A single-temperature grid cannot define a slope or bracket an inversion
      # (the ill-posed case the old LS slope returned NaN for). NA the group and
      # warn, rather than reporting a spurious single-point derivative.
      z_full <- ctmax_full <- ct1_full <- rep(NA_real_, np)
    } else {
      # ctmax_grid on the ORIGINAL temperature scale, so the inverted CTmax is
      # too (temp_mean is guaranteed non-NULL whenever CTmax / T_crit are asked).
      ctmax_grid <- if (is.null(temp_mean)) gi$tg else gi$tg + temp_mean
      zc <- tls_zct_draws(logLT[, gi$plus,  drop = FALSE],
                          logLT[, gi$minus, drop = FALSE], h, gi$tg,
                          logLT[, gi$ctmax, drop = FALSE], ctmax_grid,
                          logLT[, gi$mid0],
                          target = log_tref, target_1hr = target_1hr,
                          Tbar = temp_mean, use_closed = mid_rel_closed,
                          z_local = FALSE,
                          want_ctmax = want_ctmax, want_ct1 = want_tcrit)
      z_full <- rep(NA_real_, np)
      z_full[zc$z$draws$.draw] <- zc$z$draws$z      # pooled per-draw z (NA where non-finite)
      ctmax_full <- zc$ctmax %||% rep(NA_real_, np)
      ct1_full   <- zc$ct1   %||% rep(NA_real_, np)
    }

    if (want_z)     add("z", z_full)
    if (want_ctmax) add("CTmax", ctmax_full)
    if (want_tcrit) {
      u <- stats::runif(np, log10(TC_rate_range[1] / 100),
                        log10(TC_rate_range[2] / 100))
      add("Tcrit", ct1_full + z_full * u)
    }
  }

  if (length(allna_groups)) {
    warning("tls(): all z/CTmax draws were non-finite for group(s) ",
            paste(allna_groups, collapse = ", "),
            "; their summary rows are NA (e.g. a near-flat LT curve, a curve ",
            "that never crosses t_ref within the temperature grid, or a ",
            "single-temperature group). Do not read these as valid estimates.",
            call. = FALSE)
  }

  out <- list(
    summary = tibble::as_tibble(do.call(rbind, summ)),
    draws   = tibble::as_tibble(do.call(rbind, drw)),
    meta    = list(params = params, mode = mode, p = p, t_ref = t_ref,
                   lethal = lethal, TC_rate_range = TC_rate_range,
                   temp_mean = temp_mean, by = by)
  )
  class(out) <- c("tls", "list")
  out
}

#' @rdname tls
#' @export
tls_z <- function(object, ...) tls(object, params = "z", ...)

#' @rdname tls
#' @export
tls_ctmax <- function(object, ...) tls(object, params = "ctmax", ...)

#' @rdname tls
#' @export
tls_tcrit <- function(object, ...) tls(object, params = "tcrit", lethal = TRUE, ...)

#' @export
print.tls <- function(x, ...) {
  cat(sprintf("<tls> %s threshold; quantities: %s\n",
              x$meta$mode, paste(x$meta$params, collapse = ", ")))
  print(x$summary)
  invisible(x)
}
