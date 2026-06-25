# Accessors that expose the full posterior draws of derived TDT and
# heat-injury quantities as long-format tibbles. The motivating use case is
# group contrasts: fit one workflow per group (life stage, sex, treatment),
# run extract_tdt() / predict_*() on each, and use these helpers to pair the
# per-draw posteriors so a user can compute the posterior of (group A −
# group B) directly.

#' Posterior draws of z (thermal sensitivity, °C per 10-fold change in time)
#'
#' Pulls the per-draw values of z from an [extract_tdt()] result.
#'
#' @param et The list returned by [extract_tdt()].
#' @return A tibble with columns `.draw` and `z`.
#' @seealso [get_z_summary()] for the median + 95% credible-interval summary.
#' @examples
#' \dontrun{
#' et_a <- extract_tdt(wf_group_a)
#' et_b <- extract_tdt(wf_group_b)
#' contrast <- merge(get_z_draws(et_a), get_z_draws(et_b),
#'                   by = ".draw", suffixes = c("_a", "_b"))
#' contrast$diff <- contrast$z_a - contrast$z_b
#' quantile(contrast$diff, c(0.025, 0.5, 0.975))
#' }
#' @export
get_z_draws <- function(et) {
  stop_if_not_extract_tdt(et)
  tibble::as_tibble(et$z$draws)
}

#' Posterior summary of z (median and 95% credible interval)
#'
#' Pulls the z summary table from an [extract_tdt()] result — the headline
#' median and credible-interval numbers, as opposed to the per-draw values
#' returned by [get_z_draws()].
#'
#' @param et The list returned by [extract_tdt()].
#' @return A tibble with columns `z_median`, `z_lower`, `z_upper`.
#' @seealso [get_z_draws()] for the per-draw posterior.
#' @examples
#' \dontrun{
#' get_z_summary(extract_tdt(wf))
#' }
#' @export
get_z_summary <- function(et) {
  stop_if_not_extract_tdt(et)
  tibble::as_tibble(et$z$summary)
}

#' Posterior draws of CTmax at the reference duration
#'
#' Pulls the per-draw CTmax temperature: the temperature at which the
#' [extract_tdt()] threshold is reached after `t_ref` exposure (default 60
#' min). The threshold itself follows the `target_surv` argument of
#' [extract_tdt()] — `(low + up)/2` per draw by default, absolute 0.5 if
#' `target_surv = "absolute"`.
#'
#' @param et The list returned by [extract_tdt()].
#' @return A tibble with columns `.draw` and `CTmax` (plus the moderator
#'         column(s) for a grouped fit).
#' @seealso [get_ctmax_summary()] for the median + 95% credible-interval summary.
#' @examples
#' \dontrun{
#' get_ctmax_draws(extract_tdt(wf))
#' }
#' @export
get_ctmax_draws <- function(et) {
  stop_if_not_extract_tdt(et)
  d     <- et$CTmax$draws
  gcols <- setdiff(names(d), c(".draw", "temp"))   # moderator column(s) for a grouped fit
  out   <- tibble::tibble(.draw = d$.draw, CTmax = d$temp)
  if (length(gcols)) out <- tibble::as_tibble(cbind(d[gcols], out))
  out
}

#' Posterior summary of CTmax at the reference duration
#'
#' Pulls the CTmax summary table from an [extract_tdt()] result (the
#' temperature at which the threshold is reached after `t_ref` exposure), as
#' opposed to the per-draw values returned by [get_ctmax_draws()].
#'
#' @param et The list returned by [extract_tdt()].
#' @return A tibble with columns `temp_lower`, `temp_median`, `temp_upper`.
#' @seealso [get_ctmax_draws()] for the per-draw posterior.
#' @examples
#' \dontrun{
#' get_ctmax_summary(extract_tdt(wf))
#' }
#' @export
get_ctmax_summary <- function(et) {
  stop_if_not_extract_tdt(et)
  tibble::as_tibble(et$CTmax$summary)
}

#' Posterior draws of T_crit (rate-multiplier definition)
#'
#' Pulls the per-draw T_crit values from an [extract_tdt()] result. Errors if
#' `extract_tdt()` was called with the default `lethal = FALSE` — T_crit is
#' only meaningful for lethal-endpoint data.
#'
#' @param et The list returned by `extract_tdt(..., lethal = TRUE)`.
#' @return A tibble with columns `.draw`, `T_crit`, and `log10_rate` (the
#'         sampled log10 rate-multiplier for that posterior draw); plus the
#'         moderator column(s) for a grouped fit.
#' @seealso [get_tcrit_summary()] for the median + 95% credible-interval summary.
#' @examples
#' \dontrun{
#' get_tcrit_draws(extract_tdt(wf, lethal = TRUE))
#' }
#' @export
get_tcrit_draws <- function(et) {
  stop_if_not_extract_tdt(et)
  if (is.null(et$T_crit))
    stop("This extract_tdt() result has no T_crit. ",
         "Refit with `lethal = TRUE`.", call. = FALSE)
  d     <- et$T_crit$draws
  gcols <- setdiff(names(d), c(".draw", "temp", "log10_rate"))   # moderator(s) for a grouped fit
  out   <- tibble::tibble(.draw = d$.draw, T_crit = d$temp, log10_rate = d$log10_rate)
  if (length(gcols)) out <- tibble::as_tibble(cbind(d[gcols], out))
  out
}

#' Posterior summary of T_crit (rate-multiplier definition)
#'
#' Pulls the T_crit summary table from an [extract_tdt()] result — median and
#' credible interval, plus the rate-floor range used — as opposed to the
#' per-draw values returned by [get_tcrit_draws()]. Errors if `extract_tdt()`
#' was called with the default `lethal = FALSE`; T_crit is only meaningful for
#' lethal-endpoint data.
#'
#' @param et The list returned by `extract_tdt(..., lethal = TRUE)`.
#' @return A tibble with columns `TC_rate_low`, `TC_rate_high`, `temp_lower`,
#'         `temp_median`, `temp_upper`.
#' @seealso [get_tcrit_draws()] for the per-draw posterior.
#' @examples
#' \dontrun{
#' get_tcrit_summary(extract_tdt(wf, lethal = TRUE))
#' }
#' @export
get_tcrit_summary <- function(et) {
  stop_if_not_extract_tdt(et)
  if (is.null(et$T_crit))
    stop("This extract_tdt() result has no T_crit. ",
         "Refit with `lethal = TRUE`.", call. = FALSE)
  tibble::as_tibble(et$T_crit$summary)
}

#' Paired posterior draws of z, CTmax and (optionally) T_crit
#'
#' Merges the per-quantity draws from an [extract_tdt()] result into a single
#' tibble, joined on `.draw` (and the moderator column(s) for a grouped fit) so
#' that every quantity is on the **same posterior draw**. This preserves their
#' joint correlation, which is what per-draw contrasts, ratios and bivariate
#' plots need. `z` and `CTmax` are always returned; `T_crit` is included only
#' when [extract_tdt()] was called with `lethal = TRUE`.
#'
#' The merge is an inner join, not a column-bind: each quantity is filtered to
#' finite values independently inside [extract_tdt()] (e.g. a draw with a
#' near-zero midpoint slope has no finite `z`; a draw whose curve never crosses
#' the reference duration has no finite `CTmax`), so the draws can differ in
#' which rows they retain. Joining on `.draw` keeps only draws present for every
#' returned quantity and guarantees alignment; binding columns could silently
#' mis-pair them.
#'
#' @param et The list returned by [extract_tdt()].
#' @return A tibble with one row per shared posterior draw: columns `.draw`, the
#'   moderator column(s) for a grouped fit, `z`, `CTmax`, and (when available)
#'   `T_crit`. The auxiliary `log10_rate` from [get_tcrit_draws()] is dropped.
#' @seealso [get_z_draws()], [get_ctmax_draws()], [get_tcrit_draws()] for the
#'   individual per-quantity draws; [extract_tdt()].
#' @examples
#' \dontrun{
#' wf <- fit_4pl(std)
#' d  <- get_tls_draws(extract_tdt(wf, lethal = TRUE))
#' # joint per-draw transform keeps z and CTmax paired:
#' stats::quantile(d$CTmax - 2.5 * d$z, c(0.025, 0.5, 0.975))
#' }
#' @export
get_tls_draws <- function(et) {
  stop_if_not_extract_tdt(et)
  zd  <- get_z_draws(et)        # .draw, [moderators], z
  cd  <- get_ctmax_draws(et)    # .draw, [moderators], CTmax
  # Join keys = every shared column except the value columns (i.e. .draw plus any
  # moderator columns), so this works for single-condition and grouped fits.
  out <- dplyr::inner_join(
    zd, cd, by = intersect(setdiff(names(zd), "z"), setdiff(names(cd), "CTmax")))
  if (!is.null(et$T_crit)) {
    td <- get_tcrit_draws(et)
    td[["log10_rate"]] <- NULL  # keep .draw / moderators / T_crit
    out <- dplyr::inner_join(
      out, td, by = intersect(setdiff(names(out), c("z", "CTmax")),
                              setdiff(names(td), "T_crit")))
  }
  out
}

#' Combined posterior summary of z, CTmax and (optionally) T_crit
#'
#' Reshapes the per-quantity summaries from an [extract_tdt()] result into one
#' tidy-long tibble — the summary companion to [get_tls_draws()], and the same
#' shape [tls()] returns: one row per quantity (per group, for a grouped fit),
#' with harmonised `median` / `lower` / `upper` columns. `z` and `CTmax` are
#' always returned; `Tcrit` is included only when [extract_tdt()] was called
#' with `lethal = TRUE`. (The individual accessors use quantity-specific column
#' names — `z_median` vs `temp_median` — so this is the tidy way to tabulate or
#' plot all quantities together.)
#'
#' @param et The list returned by [extract_tdt()].
#' @return A tibble with columns `quantity` (`"z"`, `"CTmax"`, and when available
#'   `"Tcrit"`), `median`, `lower`, `upper`, plus the moderator column(s) for a
#'   grouped fit. The `T_crit` rate-floor columns (`TC_rate_*`) are dropped.
#' @seealso [get_tls_draws()] for the paired per-draw posterior;
#'   [get_z_summary()], [get_ctmax_summary()], [get_tcrit_summary()] for the
#'   individual summaries; [tls()], which returns the same tidy shape from a fit.
#' @examples
#' \dontrun{
#' get_tls_summary(extract_tdt(wf, lethal = TRUE))
#' }
#' @export
get_tls_summary <- function(et) {
  stop_if_not_extract_tdt(et)
  pick <- function(summ, quantity, med, lo, hi) {
    mods <- setdiff(names(summ),
                    c(med, lo, hi, grep("^TC_rate", names(summ), value = TRUE)))
    out  <- tibble::tibble(quantity = quantity,
                           median = summ[[med]], lower = summ[[lo]], upper = summ[[hi]])
    if (length(mods)) tibble::as_tibble(cbind(summ[mods], out)) else out
  }
  parts <- list(
    pick(get_z_summary(et),     "z",     "z_median",    "z_lower",    "z_upper"),
    pick(get_ctmax_summary(et), "CTmax", "temp_median", "temp_lower", "temp_upper"))
  if (!is.null(et$T_crit))
    parts <- c(parts, list(pick(get_tcrit_summary(et), "Tcrit",
                                "temp_median", "temp_lower", "temp_upper")))
  dplyr::bind_rows(parts)
}

#' Extract TLS estimates (draws or summary) from a `tls` object
#'
#' Pulls the posterior `draws` or the `summary` (median + 95% credible interval)
#' for any or all thermal-load-sensitivity quantities (`z`, `CTmax`, and, when the
#' fit was extracted with `lethal = TRUE`, `Tcrit`) from a [tls()] result. This is
#' the accessor companion for `tls` objects, mirroring [get_tls_summary()] /
#' [get_tls_draws()] for an [extract_tdt()] result. Draws are returned in
#' tidy-long form (`quantity`, `.draw`, `value`, plus any grouping column), so a
#' posterior contrast is a join on `.draw` followed by a difference of `value`.
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
#' @seealso [tls()]; [get_tls_summary()] and [get_tls_draws()] for an
#'   [extract_tdt()] result.
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
    stop("`x` must be a `tls` object (the result of `tls()`). For an ",
         "`extract_tdt()` result use `get_tls_summary()` / `get_tls_draws()`.",
         call. = FALSE)
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

# Internal: typeguard for extract_tdt() results.
stop_if_not_extract_tdt <- function(x) {
  if (!is.list(x) || is.null(x$z) || is.null(x$CTmax))
    stop("Expected an extract_tdt() result.", call. = FALSE)
  invisible(NULL)
}
