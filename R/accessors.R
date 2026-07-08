# Accessors that expose the full posterior draws of derived TLS and
# heat-injury quantities as long-format tibbles. The motivating use case is
# group contrasts: fit one workflow per group (life stage, sex, treatment),
# run tls() / predict_*() on each, and use these helpers to pair the
# per-draw posteriors so a user can compute the posterior of (group A −
# group B) directly.

#' Extract TLS estimates (draws or summary) from a `tls` object
#'
#' Pulls the posterior `draws` or the `summary` (median + 95% credible interval)
#' for any or all thermal-load-sensitivity quantities (`z`, `CTmax`, and, when the
#' fit was extracted with `lethal = TRUE`, `Tcrit`) from a [tls()] result. Draws
#' are returned in tidy-long form (`quantity`, `.draw`, `value`, plus any grouping
#' column), so a posterior contrast is a join on `.draw` followed by a difference
#' of `value`.
#'
#' @param x A `tls` object, the result of [tls()].
#' @param what Either `"summary"` (the default; median + lower/upper per
#'   quantity x group) or `"draws"` (one row per posterior draw per quantity x
#'   group).
#' @param params Optional character vector selecting quantities to return, matched
#'   case-insensitively against `"z"`, `"CTmax"`, `"Tcrit"`. `NULL` (default)
#'   returns all available quantities.
#' @return A tibble. For `what = "summary"`: columns `quantity`, `median`,
#'   `lower`, `upper` (plus grouping columns for a grouped fit). For
#'   `what = "draws"`: columns `quantity`, `.draw`, `value` (plus grouping
#'   columns).
#' @seealso [tls()], which produces the object this reads.
#' @examples
#' \dontrun{
#' fit <- fit_4pl(standardize_data(my_data, temp = "temp", duration = "time",
#'                                 n_total = "n", n_dead = "dead"))
#' est <- tls(fit, lethal = TRUE)
#' get_tls_est(est, "summary")             # all quantities: median + 95% CrI
#' zd  <- get_tls_est(est, "draws", "z")   # z draws, ready for a contrast
#' }
#' @export
get_tls_est <- function(x, what = c("summary", "draws"), params = NULL) {
  if (!inherits(x, "tls"))
    stop("`x` must be a `tls` object (the result of `tls()`).", call. = FALSE)
  what <- match.arg(what)
  out  <- if (what == "summary") x$summary else x$draws
  if (!is.null(params)) {
    avail <- unique(out$quantity)
    sel   <- avail[tolower(avail) %in% tolower(params)]
    if (!length(sel))
      stop("No matching `params`. Available quantities: ",
           paste(avail, collapse = ", "), ".", call. = FALSE)
    out <- out[out$quantity %in% sel, , drop = FALSE]
  }
  tibble::as_tibble(out)
}

#' Posterior estimates of the natural-scale 4PL parameters (draws or summary)
#'
#' The 4PL-shape companion to [get_tls_est()]. Returns the posterior of the
#' natural-scale 4PL curve parameters from a fitted [fit_4pl()] workflow, as
#' either per-draw values or a median + 95% credible-interval summary, per
#' moderator group. Whereas [get_tls_est()] operates on a [tls()] object and
#' returns the derived TLS quantities (z, CTmax, Tcrit), this operates on the
#' fitted workflow and returns the raw curve parameters those quantities are
#' derived from. A thin wrapper: `"summary"` calls [tdt_parameter_table()],
#' `"draws"` calls [extract_4pl_pars()].
#'
#' @param workflow A fitted [fit_4pl()] workflow (`workflow$fit` not `NULL`).
#' @param what `"summary"` (default; one row per parameter with median + 95% CrI
#'   at the centring temperature `T_bar`, via [tdt_parameter_table()] -- this
#'   also includes the derived `z` slope row) or `"draws"` (per-draw natural-
#'   scale `low`, `up`, `k`, `mid` at each assay temperature, via
#'   [extract_4pl_pars()]).
#' @param by Optional moderator column(s) to group by, passed through. `NULL`
#'   (default) uses the fit's own grouping for a grouped fit, else one group.
#' @param temps For `what = "draws"`, numeric vector of assay temperatures (°C)
#'   to evaluate the curve parameters at, passed to [extract_4pl_pars()].
#'   `NULL` (default) uses the fit's observed unique assay temperatures.
#'   Ignored for `what = "summary"` (which always reports the `T_bar` slice).
#' @return A tibble (plus the moderator column(s) for a grouped fit): for
#'   `"summary"`, one row per parameter with median and credible bounds at
#'   `T_bar`; for `"draws"`, one row per posterior draw × temperature with
#'   `temp`, `low`, `up`, `k`, `mid` (and `.draw`).
#' @seealso [get_tls_est()] for the derived TLS quantities; [tdt_parameter_table()]
#'   and [extract_4pl_pars()], which this wraps.
#' @examples
#' \dontrun{
#' wf <- fit_4pl(standardize_data(my_data, temp = "temp", duration = "time",
#'                                n_total = "n", n_dead = "dead"))
#' get_4pl_est(wf, "summary")                 # low/up/k/mid median + 95% CrI
#' get_4pl_est(wf, "draws", by = "species")   # per-draw params, per group
#' }
#' @export
get_4pl_est <- function(workflow, what = c("summary", "draws"), by = NULL,
                        temps = NULL) {
  what <- match.arg(what)
  if (what == "summary") tdt_parameter_table(workflow, by = by)
  else extract_4pl_pars(workflow, temps = temps, by = by)
}

#' Posterior draws of heat-injury trajectory
#'
#' Long-format tibble of per-draw HI, dose, survival and mortality at every
#' time step of the supplied temperature trace. Requires that
#' [predict_heat_injury()] was called with `save_draws = TRUE`.
#'
#' @param hi The list returned by `predict_heat_injury(..., save_draws = TRUE)`.
#' @return A tibble with columns `.draw`, `time`, `temp`, `dose`, `hi`,
#'         `survival`, `mortality`.
#' @examples
#' \dontrun{
#' hi <- predict_heat_injury(trace, wf, save_draws = TRUE)
#' get_hi_draws(hi)
#' }
#' @export
get_hi_draws <- function(hi) {
  if (!is.list(hi) || is.null(hi$summary) ||
      !"time" %in% names(hi$summary))
    stop("Expected a predict_heat_injury() result.", call. = FALSE)
  if (is.null(hi$draws))
    stop("hi$draws is NULL. ",
         "Re-run predict_heat_injury() with `save_draws = TRUE`.",
         call. = FALSE)
  tibble::as_tibble(hi$draws)
}

#' Posterior draws of survival
#'
#' Accepts either a [predict_survival_curves()] result (static survival on a
#' temperature × duration grid) or a [predict_heat_injury()] result (dynamic
#' survival along a temperature trace). In the heat-injury case requires
#' `save_draws = TRUE`.
#'
#' @param x A list returned by [predict_survival_curves()] or by
#'          [predict_heat_injury()].
#' @return A long-format tibble. For [predict_survival_curves()]: `.draw`,
#'         `temp`, `duration`, `survival`. For [predict_heat_injury()]:
#'         `.draw`, `time`, `temp`, `survival`. For a grouped fit the moderator
#'         column(s) are also carried (so per-group draws are not collapsed).
#' @examples
#' \dontrun{
#' psc <- predict_survival_curves(wf, temps = c(32, 34),
#'                                durations = c(0.5, 1))
#' get_surv_draws(psc)
#'
#' hi <- predict_heat_injury(trace, wf, save_draws = TRUE)
#' get_surv_draws(hi)
#' }
#' @export
get_surv_draws <- function(x) {
  if (is.list(x) && !is.null(x$draws_matrix) && !is.null(x$grid)) {
    return(surv_grid_to_long(x))
  }
  if (is.list(x) && !is.null(x$summary) &&
      "time" %in% names(x$summary)) {
    if (is.null(x$draws))
      stop("x$draws is NULL. ",
           "Re-run predict_heat_injury() with `save_draws = TRUE`.",
           call. = FALSE)
    # Carry any moderator column(s) for a grouped fit so per-group survival draws
    # are not collapsed onto colliding (.draw, time) keys. integrate_pars()
    # prepends the `by` column(s); anything beyond the standard heat-injury draw
    # columns is a moderator. For a single-condition fit `gcols` is empty and the
    # output is unchanged.
    gcols <- setdiff(names(x$draws),
                     c(".draw", "time", "temp", "dose", "hi", "survival", "mortality"))
    return(tibble::as_tibble(
      x$draws[, c(gcols, ".draw", "time", "temp", "survival"), drop = FALSE]
    ))
  }
  stop("Could not recognize input as a predict_survival_curves() or ",
       "predict_heat_injury() result.", call. = FALSE)
}

# Internal: predict_survival_curves draws_matrix [ndraws x ngrid] → long. Carries
# the moderator column(s) (those in $summary beyond temp/duration/survival_*) for
# a grouped fit, so per-group draws are not collapsed onto colliding (temp,
# duration) keys; for a single-condition fit there are none and the output is
# unchanged.
surv_grid_to_long <- function(psc) {
  mat  <- psc$draws_matrix
  grid <- psc$grid
  if (ncol(mat) != nrow(grid))
    stop("draws_matrix and grid dimensions disagree.", call. = FALSE)
  gcols <- setdiff(names(psc$summary),
                   c("temp", "duration", "survival_lower", "survival_median", "survival_upper"))
  ndraws <- nrow(mat)
  tibble::as_tibble(do.call(rbind, lapply(seq_len(ncol(mat)), function(j) {
    base <- data.frame(.draw    = seq_len(ndraws),
                       temp     = grid$temp[j],
                       duration = grid$duration[j],
                       survival = mat[, j])
    if (length(gcols)) base <- cbind(grid[j, gcols, drop = FALSE], base, row.names = NULL)
    base
  })))
}
