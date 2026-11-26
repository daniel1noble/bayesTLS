# Generic, moderator-aware TDT extraction for arbitrary fitted 4PL models.
# tls() evaluates the four 4PL sub-parameters at a user-defined moderator x
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
#'   midpoint. The dose-response model is read from the workflow; a bare
#'   `brmsfit` is recognised from its formula, and a hand-written formula is
#'   treated as the log-logistic.
#' @param by Character vector of moderator columns defining the groups reported
#'   separately (e.g. `"sex"`). `NULL` (default) uses the fit's own groups, as
#'   the other readers do: a `bayes_tls` workflow fitted with a moderator
#'   (e.g. `ctmax = ~ 0 + recovery`, or `by =` in [fit_4pl()]) is reported per
#'   group, and a single-condition fit as one group. A bare `brmsfit` whose
#'   CTmax/z (or midpoint) vary by a moderator needs `by` given explicitly
#'   (an error says so). With a user-supplied `newdata` grid, `NULL` treats
#'   the whole grid as one group.
#' @param params Quantities to return: `"all"` (z, CTmax, and T_crit when
#'   `lethal = TRUE`), or a subset of `c("z", "ctmax", "tcrit")`.
#' @param target_surv Survival threshold the LT curve is read at:
#'   `"relative"` (default; the per-draw midpoint
#'   `(low + up)/2`), `"absolute"` (literal LT50), or a numeric `p` in `(0, 1)`
#'   for an absolute LT`p`.
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
#' @param temp_grid Temperatures **on the centred `temp` scale** (temperature
#'   minus `temp_mean`, unlike the derive_*() readers, which take degrees C) at
#'   which the LT curve is evaluated. `z` is the local `z(T)` pooled over this
#'   grid. `CTmax` (and the 1-hour anchor of `T_crit`) is exact when the LT
#'   curve is a straight line in temperature; when it bends (see Details) it is
#'   inverted on this grid by linear interpolation between its points (a finer
#'   grid reduces that error), and draws that do not cross on it are searched
#'   outward, beyond each end by the data's temperature range or the grid's own
#'   width, whichever is larger. Draws left without a crossing are reported in
#'   a warning. Default: 11 points over the observed range. A grid lying
#'   entirely outside the observed temperatures, or reaching more than the
#'   data's own temperature range beyond them, warns, since that is what a
#'   grid given in degrees C usually looks like (the check cannot tell the
#'   scales apart for data near 0 degrees C).
#' @param re_formula Passed to [brms::posterior_linpred()]; default `NA`
#'   (fixed effects only) -- EXCEPT when a `by` variable enters the model only
#'   as a random-effect grouping factor (e.g. a phylogenetic
#'   `(1 | gr(species, cov = A))` term): the default then auto-includes
#'   exactly the matching random-effect term(s), with a message, so
#'   per-level estimates do not silently collapse to the population mean.
#'   Passing `re_formula` explicitly (even `NA`) always overrides the
#'   auto-rule.
#'   (population level). Include group-level terms (e.g. clone) for per-group draws.
#' @param lower,upper Response bounds of the disjoint-bounds reparameterisation.
#'   Read from a `bayes_tls` workflow (the bounds it was fitted with); supply
#'   them for a bare `brmsfit` fitted with non-default bounds, or to override.
#' @param nlpars Non-linear parameter names, ordered low / up / k / mid (raw scale).
#' @param ndraws Optional number of posterior draws to subsample: one
#'   subsample, shared by all four sub-parameters. A value at least the
#'   posterior size uses all draws.
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
#' `tls()` derives `z`, `CTmax` and `T_crit` from the same posterior draws, so
#' on a common temperature grid they are mutually consistent for every
#' threshold. `z` is the central-difference local `z(T)` (negative reciprocal of
#' the local slope of `log10(LT)` on temperature) pooled over `temp_grid`;
#' `CTmax` is the per-draw temperature at which the fitted LT curve crosses
#' `t_ref`: from the exact closed form whenever `log10(LT)` is a straight line
#' in temperature (checked on the evaluated curve, for every draw), and by
#' per-draw inversion of the bent curve otherwise. Under
#' `target_surv = "relative"` (the default) `log10(LT) = mid(T)` is linear for
#' midpoint and direct fits alike, so `z` is constant in temperature and
#' equals `-1 / b_mid_temp_c`, and `CTmax` is exact even far outside the data;
#' an absolute threshold is linear too when the asymptotes and steepness do
#' not change with temperature. With temperature effects on the asymptotes or
#' steepness, an absolute threshold bends the curve: the local `z(T)` follows
#' it, and `CTmax` is the curve's crossing of `t_ref`, interpolated linearly
#' between grid points. The crossing taken is the one where the LT shortens
#' through `t_ref` as temperature rises; a draw with none, or more than one,
#' has no `CTmax`. Near an asymptote the LT at an absolute `p` can be
#' undefined at some temperatures (`p` outside the fitted asymptotes there);
#' it is then treated as zero (upper asymptote at or below `p`) or infinite
#' (lower asymptote at or above `p`), so a crossing next to such a
#' temperature is still found.
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
  # `mode`/`p` are superseded by `target_surv`; map them when supplied so
  # existing calls keep working. target_surv defaults to "relative" (the
  # per-draw midpoint threshold).
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
  # Asymptote windows of the disjoint-bounds reparameterisation. A workflow
  # carries the ones it was fitted with (meta$bounds), as derive_z(),
  # extract_4pl_pars(), hi_pars() and tdt_parameter_table() already use;
  # explicit lower/upper override them (a missing one is filled from the
  # workflow); a bare brmsfit uses the arguments. Only absolute thresholds
  # depend on this (the relative LT is the midpoint).
  bounds <- if (missing(lower) && missing(upper) && !is.null(meta$bounds)) {
    meta$bounds
  } else {
    compute_4pl_bounds(if (missing(lower)) meta$lower %||% lower else lower,
                       if (missing(upper)) meta$upper %||% upper else upper)
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

  # `by = NULL` means the fit's own groups (it used to take the moderator from
  # data row 1 and return that one group, unlabelled).
  by <- tls_resolve_by(object, fit, by, newdata)

  mdata <- fit$data
  if (!temp %in% names(mdata)) {
    stop("Temperature column `", temp, "` not found in the model data.",
         call. = FALSE)
  }
  temp_grid <- tdt_finite_grid(temp_grid)
  tls_check_temp_grid(temp_grid, mdata[[temp]], temp_mean)

  # Auto-resolve re_formula when a `by` variable enters the model ONLY as a
  # random-effect grouping factor (e.g. phylogenetic fits with
  # ctmax = ~ 1 + (1 | gr(species, cov = A)); issue #12). With the legacy
  # hard-coded NA, posterior_linpred would exclude those effects and every
  # `by` level would silently collapse to the population mean. Only the
  # DEFAULT is resolved -- an explicit re_formula (even NA) always wins -- and
  # only the matching bar terms are included, never a blanket NULL (the grid
  # pads unrelated grouping columns with arbitrary first values, which blanket
  # inclusion would turn into phantom effects).
  if (missing(re_formula) && !is.null(by) &&
      !is.null(fit$formula$pforms)) {
    pforms_rhs <- vapply(fit$formula$pforms, function(f)
      paste(deparse(f[[length(f)]]), collapse = " "), character(1))
    fixed_vars <- unique(unlist(lapply(pforms_rhs, rhs_fixed_vars)))
    re_auto    <- tdt_auto_re_formula(pforms_rhs, by, fixed_vars)
    if (inherits(re_auto, "formula")) {
      re_formula <- re_auto
      message("`by = ", paste(by, collapse = ", "), "` enters the model only ",
              "as a random effect; including ",
              paste(deparse(re_auto), collapse = " "),
              " in posterior_linpred. Pass `re_formula` explicitly to ",
              "override.")
    }
  }

  # --- build the moderator x temperature grid (shared engine) ----------------
  base_grid <- tls_build_grid(mdata, by = by, temp = temp,
                              temp_grid = temp_grid, newdata = newdata,
                              censored = tdt_fit_censored(fit))

  # tls() derives z / CTmax / T_crit through the shared tls_zct_draws() engine:
  # the central-difference local z pooled over the grid, and either the exact
  # closed form (when logLT is a straight line in temperature) or a true
  # per-draw inversion of the bent LT curve. To do that we evaluate logLT at,
  # per group, temp_c = 0 (the closed-form anchor), the grid +/- h (local z),
  # and the grid itself (the CTmax inversion). A bent curve's draws that do not
  # cross t_ref on the grid get a second, outward search beyond it
  # (tls_extend_crossings()). This replaces the earlier single global
  # least-squares slope, which was exact only for a linear (relative) LT curve
  # and a ~linear approximation of a bent (absolute) one.
  h          <- 1e-3
  log_tref   <- log10(t_ref / time_multiplier)        # CTmax target (model units)
  target_1hr <- log10(60   / time_multiplier)         # T_crit 1-hour anchor
  # The outward search reaches the data's temperature range (or the grid's own
  # width, if wider) beyond each end of the grid.
  obs_width  <- diff(range(mdata[[temp]], na.rm = TRUE))

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
      tg = tg, L = L, modrow = modrow,
      mod   = if (is.null(by)) NULL else base_grid[rows[1], by, drop = FALSE],
      mid0  = offset + 1L,
      minus = offset + 1L + seq_len(L),
      plus  = offset + 1L + L + seq_len(L),
      ctmax = offset + 1L + 2L * L + seq_len(L))
    offset <- offset + 1L + 3L * L
  }
  aug <- do.call(rbind, aug_parts)
  rownames(aug) <- NULL

  # Reproducibility: seed the draw subsample (tls_draw_ids(), when `ndraws` is
  # set) and the T_crit rate draws below from one stream. The augmented grid adds
  # rows, not draws, so the subsample -- and hence seeded T_crit -- is unchanged.
  if (!is.null(seed)) set.seed(seed)
  # One draw subsample shared by the four sub-parameters. Passing `ndraws`
  # straight to posterior_linpred() drew a fresh random subsample for each of
  # the four (low, up, k, mid), so they came from mismatched draws (which
  # mattered for absolute thresholds).
  draw_ids <- tls_draw_ids(fit, ndraws)

  # --- evaluate sub-parameters at every grid row (shared engine) -------------
  sp    <- tls_eval_subpars(fit, aug, bounds,
                            nlpars = nlpars, re_formula = re_formula,
                            draw_ids = draw_ids, mode = mode, p = p)
  logLT <- sp$logLT
  np    <- nrow(logLT)
  p_abs <- if (identical(mode, "absolute")) p else NA_real_
  model <- tdt_model(fit)

  want_z     <- "z"     %in% params
  want_ctmax <- "ctmax" %in% params
  want_tcrit <- "tcrit" %in% params

  # --- per-group z / CTmax / 1-hour anchor from the shared draws -------------
  res  <- list()   # per group: z, ctmax, ct1, closed
  bent <- list()   # groups whose bent curve needs the outward second pass
  for (g in groups) {
    gi <- grp_cols[[g]]
    if (gi$L < 2L) {
      # A single-temperature grid cannot define a slope or bracket an inversion
      # (the ill-posed case the old LS slope returned NaN for). NA the group and
      # warn, rather than reporting a spurious single-point derivative.
      res[[g]] <- list(z = rep(NA_real_, np), ctmax = rep(NA_real_, np),
                       ct1 = rep(NA_real_, np), closed = NA)
      next
    }
    # ctmax_grid on the ORIGINAL temperature scale, so the inverted CTmax is
    # too (temp_mean is guaranteed non-NULL whenever CTmax / T_crit are asked).
    ctmax_grid <- if (is.null(temp_mean)) gi$tg else gi$tg + temp_mean
    spc <- lapply(sp[c("low", "up", "k", "mid", "logLT")],
                  function(m) m[, gi$ctmax, drop = FALSE])
    zc <- tls_zct_draws(logLT[, gi$plus,  drop = FALSE],
                        logLT[, gi$minus, drop = FALSE], h, gi$tg,
                        spc$logLT, ctmax_grid, logLT[, gi$mid0],
                        target = log_tref, target_1hr = target_1hr,
                        Tbar = temp_mean, z_local = FALSE,
                        want_ctmax = want_ctmax, want_ct1 = want_tcrit,
                        sp_ctmax = spc, p = p_abs, model = model)
    z_full <- rep(NA_real_, np)
    z_full[zc$z$draws$.draw] <- zc$z$draws$z      # pooled per-draw z (NA where non-finite)
    res[[g]] <- list(z = z_full, ctmax = zc$ctmax %||% rep(NA_real_, np),
                     ct1 = zc$ct1 %||% rep(NA_real_, np), closed = zc$closed)
    # A bent curve was inverted on the grid only: the draws that do not cross
    # there get an outward search beyond it, all groups together.
    if (isFALSE(zc$closed)) {
      keep <- c(if (want_ctmax) "ctmax", if (want_tcrit) "ct1")
      bent[[g]] <- list(newrow = gi$modrow, tg = gi$tg, sp = spc,
                        crossings = res[[g]][keep],
                        targets = list(ctmax = log_tref, ct1 = target_1hr)[keep])
    }
  }
  if (length(bent)) {
    bent <- tls_extend_crossings(fit, bent, temp, width = obs_width,
                                 Tbar = temp_mean %||% 0, draw_ids = draw_ids,
                                 bounds = bounds, nlpars = nlpars,
                                 re_formula = re_formula, mode = mode, p = p_abs,
                                 model = model)
    for (g in names(bent)) res[[g]][names(bent[[g]]$crossings)] <- bent[[g]]$crossings
  }

  # --- summarise per group (T_crit rate draws in group order, as before) -----
  summ <- list(); drw <- list()
  allna_groups <- character(0)   # groups whose z/CTmax draws are all non-finite
  uncrossed    <- character(0)   # groups with draws whose crossing was not found
  for (g in groups) {
    gi    <- grp_cols[[g]]
    gcols <- gi$mod
    add <- function(q, v) {
      # Summarise on the finite draws only. A near-flat LT curve sends
      # z = -1/slope to +/-Inf, a curve with no crossing of its target gives an
      # NA CTmax, and a single-temperature group is NA'd above; without this,
      # stats::quantile() aborts ("missing values not allowed") and takes down
      # the whole tls() call. The raw draws (Inf/NaN included) are still
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
    r <- res[[g]]
    if (want_z)     add("z", r$z)
    if (want_ctmax) add("CTmax", r$ctmax)
    if (want_tcrit) {
      u <- stats::runif(np, log10(TC_rate_range[1] / 100),
                        log10(TC_rate_range[2] / 100))
      add("Tcrit", r$ct1 + r$z * u)
    }
    # Report the draws a bent curve left without a crossing, per quantity and
    # target (a quantity with no finite draw at all is reported below instead).
    if (isFALSE(r$closed)) {
      miss <- function(v) { n <- sum(!is.finite(v)); if (n < np) n else 0L }
      parts <- c(if (want_ctmax && miss(r$ctmax) > 0L)
                   sprintf("CTmax %d of %d draws (t_ref = %g min)", miss(r$ctmax), np, t_ref),
                 if (want_tcrit && miss(r$ct1) > 0L)
                   sprintf("T_crit %d of %d draws (its 1-hour anchor)", miss(r$ct1), np))
      if (length(parts)) {
        ext <- range(tls_extend_grid(gi$tg, width = obs_width)) + (temp_mean %||% 0)
        uncrossed <- c(uncrossed, sprintf("%s: %s; searched %.4g to %.4g",
                                          g, paste(parts, collapse = ", "),
                                          ext[1], ext[2]))
      }
    }
  }

  if (length(uncrossed)) {
    warning("tls(): some draws have no single crossing of their target exposure ",
            "within the temperatures searched, so they are NA and left out of ",
            "the summary, which can then under-represent the extremes. ",
            paste(uncrossed, collapse = "; "), ". A draw has none when its ",
            "curve stays on one side of the target over the search (the crossing ",
            "lies further out, or the target is outside the fitted asymptotes ",
            "there) or crosses it more than once. The search covers temp_grid ",
            "extended beyond each end by the data's temperature range or the ",
            "grid's own width, whichever is larger; a wider temp_grid (centred ",
            "scale) helps only when the crossing lies further out.",
            call. = FALSE)
  }
  if (length(allna_groups)) {
    warning("tls(): all z/CTmax draws were non-finite for group(s) ",
            paste(allna_groups, collapse = ", "),
            "; their summary rows are NA (e.g. a near-flat LT curve, a curve ",
            "with no crossing of its target within the temperatures searched, ",
            "or a single-temperature group). Do not read these as valid ",
            "estimates.", call. = FALSE)
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

# Resolve tls()'s `by`: an explicit `by` wins; a user-supplied `newdata` grid is
# taken as given (one group unless `by` is set); a bayes_tls workflow uses its
# recorded moderators via tdt_resolve_by() (as derive_z() and the other readers
# do), so a grouped fit is reported per group; a bare brmsfit that looks grouped
# (multi-level CTmax/z or mid coefficients) must name its moderator(s), since
# filling the moderator from data row 1 would silently report one group.
tls_resolve_by <- function(object, fit, by, newdata = NULL) {
  if (!is.null(by) || !is.null(newdata)) return(by)
  if (inherits(object, "bayes_tls")) return(tdt_resolve_by(object, by))
  wrapped <- structure(list(fit = fit, meta = list()), class = "bayes_tls")
  if (tdt_is_grouped(wrapped))
    stop("This brmsfit appears to vary CTmax/z (or mid) by a moderator. Pass ",
         "`by = <column(s)>` to tls() to extract per group.", call. = FALSE)
  NULL
}

# Warn when a tls() temp_grid looks like degrees C rather than centred
# temperatures: it lies entirely outside the observed (centred) temperatures, or
# reaches more than the data's own temperature range beyond them. tls() reads
# temp_grid on the centred scale, so a grid given in degrees C (as the
# derive_*() readers take it) lands off the data and silently returns nonsense
# z / CTmax. A centred grid that overlaps the data and extends moderately
# beyond it (e.g. to bracket an absolute CTmax) stays silent. The two scales
# cannot be told apart when the centring temperature is small next to the
# data's temperature range (data near 0 degrees C), so the check can miss
# those.
tls_check_temp_grid <- function(temp_grid, temp_obs, temp_mean = NULL) {
  if (is.null(temp_grid) || !length(temp_grid)) return(invisible(FALSE))
  rng   <- range(temp_obs, na.rm = TRUE)
  w     <- max(diff(rng), 1)
  gr    <- range(temp_grid, na.rm = TRUE)
  apart <- gr[1] > rng[2] || gr[2] < rng[1]
  far   <- gr[1] < rng[1] - w || gr[2] > rng[2] + w
  if (!isTRUE(apart) && !isTRUE(far)) return(invisible(FALSE))
  hint <- if (!is.null(temp_mean) && is.finite(temp_mean))
    sprintf(" If you meant degrees C, pass temp_grid - %.4g (the centring temperature).",
            temp_mean) else ""
  warning(sprintf(paste0(
    "`temp_grid` (%.4g to %.4g) %s the observed temperatures (%.4g to %.4g on ",
    "the centred scale). tls() reads temp_grid as centred temperatures ",
    "(temperature minus the centring temperature), so the z and CTmax read on ",
    "this grid are extrapolations.%s"),
    gr[1], gr[2], if (isTRUE(apart)) "lies entirely outside" else
      "reaches far outside", rng[1], rng[2], hint), call. = FALSE)
  invisible(TRUE)
}
