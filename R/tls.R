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
#' @param mode Threshold for the LT curve: `"absolute"` (default; the `p`
#'   survival LT) or `"relative"` (midpoint between the fitted asymptotes).
#' @param p Survival threshold for `mode = "absolute"`. Default 0.5.
#' @param t_ref Reference exposure duration for `CTmax`, in minutes. Default 60.
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
#'   evaluated for the slope/crossing. Default: 11 points over the observed range.
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
#' @examples
#' \dontrun{
#' tls(joint_sex_fit, by = "sex", lethal = TRUE, temp_mean = 36.1)  # z, CTmax, T_crit per sex
#' tls(wf_leaf, params = "z")                                       # z only, workflow
#' }
#' @export
tls <- function(object, by = NULL, params = "all",
                mode = c("absolute", "relative"), p = 0.5,
                t_ref = 60, time_multiplier = NULL,
                lethal = FALSE, TC_rate_range = c(0.1, 1),
                temp = "temp_c", temp_mean = NULL, temp_grid = NULL,
                re_formula = NA, lower = 0, upper = 1,
                nlpars = c("lowraw", "upraw", "logk", "mid"),
                ndraws = NULL, newdata = NULL,
                probs = c(0.025, 0.5, 0.975), seed = NULL) {
  mode <- match.arg(mode)

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

  # --- build the moderator x temperature grid --------------------------------
  if (is.null(newdata)) {
    if (is.null(temp_grid)) {
      temp_grid <- seq(min(mdata[[temp]]), max(mdata[[temp]]), length.out = 11)
    }
    base <- if (is.null(by)) data.frame(.dummy = 1) else unique(mdata[, by, drop = FALSE])
    newdata <- do.call(rbind, lapply(seq_len(nrow(base)), function(i) {
      g <- base[rep(i, length(temp_grid)), , drop = FALSE]
      g[[temp]] <- temp_grid
      g
    }))
    newdata$.dummy <- NULL
    rownames(newdata) <- NULL
  }
  # Fill any model-data columns the grid lacks so posterior_linpred validates.
  for (v in setdiff(names(mdata), names(newdata))) {
    newdata[[v]] <- mdata[[v]][1]
  }
  newdata$.grp <- if (is.null(by)) "all" else
    do.call(paste, c(newdata[by], sep = " / "))

  # Reproducibility: seed the posterior_linpred draw subsample (when `ndraws` is
  # set) and the T_crit rate draws below from one stream.
  if (!is.null(seed)) set.seed(seed)

  # --- evaluate sub-parameters at every grid row -----------------------------
  b  <- compute_4pl_bounds(lower, upper)
  lp <- function(np) brms::posterior_linpred(fit, newdata = newdata, nlpar = np,
                                             re_formula = re_formula, ndraws = ndraws)
  low <- b$low_min + stats::plogis(lp(nlpars[1])) * b$low_w
  up  <- b$up_min  + stats::plogis(lp(nlpars[2])) * b$up_w
  k   <- exp(lp(nlpars[3]))
  mid <- lp(nlpars[4])
  logLT <- if (mode == "absolute") mid + log((up - p) / (p - low)) / k else mid

  # --- derive per-group quantities from the shared draws ---------------------
  tc <- newdata[[temp]]
  log_tref <- log10(t_ref / time_multiplier)
  summ <- list(); drw <- list()
  for (g in unique(newdata$.grp)) {
    cols <- which(newdata$.grp == g)
    w    <- tc[cols] - mean(tc[cols])
    slope <- as.vector(logLT[, cols, drop = FALSE] %*% w / sum(w^2))   # per-draw LS slope
    inter <- rowMeans(logLT[, cols, drop = FALSE]) - slope * mean(tc[cols])  # at temp = 0
    gcols <- if (is.null(by)) NULL else newdata[cols[1], by, drop = FALSE]
    add <- function(q, v) {
      s <- data.frame(quantity = q, median = stats::median(v),
                      lower = stats::quantile(v, probs[1], names = FALSE),
                      upper = stats::quantile(v, probs[3], names = FALSE),
                      row.names = NULL)
      d <- data.frame(quantity = q, .draw = seq_along(v), value = v,
                      row.names = NULL)
      if (!is.null(gcols)) {                       # data.frame recycles the 1-row gcols
        s <- data.frame(gcols, s, row.names = NULL)
        d <- data.frame(gcols, d, row.names = NULL)
      }
      summ[[length(summ) + 1L]] <<- s
      drw[[length(drw) + 1L]]   <<- d
    }
    z <- -1 / slope
    if ("z"     %in% params) add("z", z)
    if ("ctmax" %in% params) add("CTmax", temp_mean + (log_tref - inter) / slope)
    if ("tcrit" %in% params) {
      ct1 <- temp_mean + (log10(60 / time_multiplier) - inter) / slope  # CTmax at 1 h
      u   <- stats::runif(length(z), log10(TC_rate_range[1] / 100),
                          log10(TC_rate_range[2] / 100))
      add("Tcrit", ct1 + z * u)
    }
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
