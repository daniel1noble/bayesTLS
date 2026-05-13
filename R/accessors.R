# Accessors that expose the full posterior draws of derived TDT and
# heat-injury quantities as long-format tibbles. The motivating use case is
# group contrasts: fit one workflow per group (life stage, sex, treatment),
# run extract_tdt() / predict_*() on each, and use these helpers to pair the
# per-draw posteriors so a user can compute the posterior of (group A −
# group B) directly.

#' Posterior draws of z (TDT slope, °C per log10-time)
#'
#' Pulls the per-draw values of z from an [extract_tdt()] result.
#'
#' @param et The list returned by [extract_tdt()].
#' @return A tibble with columns `.draw` and `z`.
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

#' Posterior draws of CTmax at the reference duration
#'
#' Pulls the per-draw CTmax temperature (the temperature at which
#' `target_surv = 0.5` is reached at the [extract_tdt()] `t_ref`, default 60
#' min).
#'
#' @param et The list returned by [extract_tdt()].
#' @return A tibble with columns `.draw` and `CTmax`.
#' @examples
#' \dontrun{
#' get_ctmax_draws(extract_tdt(wf))
#' }
#' @export
get_ctmax_draws <- function(et) {
  stop_if_not_extract_tdt(et)
  d <- et$CTmax$draws
  tibble::tibble(.draw = d$.draw, CTmax = d$temp)
}

#' Posterior draws of T_crit (rate-multiplier definition)
#'
#' Pulls the per-draw T_crit values from an [extract_tdt()] result. Errors if
#' `extract_tdt()` was called with the default `lethal = FALSE` — T_crit is
#' only meaningful for lethal-endpoint data.
#'
#' @param et The list returned by `extract_tdt(..., lethal = TRUE)`.
#' @return A tibble with columns `.draw`, `T_crit`, and `log10_rate` (the
#'         sampled log10 rate-multiplier for that posterior draw).
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
  d <- et$T_crit$draws
  tibble::tibble(.draw = d$.draw, T_crit = d$temp, log10_rate = d$log10_rate)
}

#' Posterior draws of heat-injury trajectory
#'
#' Long-format tibble of per-draw HI, dose, survival and mortality at every
#' time step of the supplied temperature trace. Requires that
#' [predict_heat_injury()] was called with `save_draws = TRUE`.
#'
#' @param hi The list returned by `predict_heat_injury(..., save_draws = TRUE)`.
#' @return A tibble with columns `.draw`, `time_h`, `temp`, `dose`, `hi`,
#'         `survival`, `mortality`.
#' @examples
#' \dontrun{
#' hi <- predict_heat_injury(trace, wf, save_draws = TRUE)
#' get_hi_draws(hi)
#' }
#' @export
get_hi_draws <- function(hi) {
  if (!is.list(hi) || is.null(hi$summary) ||
      !"time_h" %in% names(hi$summary))
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
#'         `.draw`, `time_h`, `temp`, `survival`.
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
      "time_h" %in% names(x$summary)) {
    if (is.null(x$draws))
      stop("x$draws is NULL. ",
           "Re-run predict_heat_injury() with `save_draws = TRUE`.",
           call. = FALSE)
    return(tibble::as_tibble(
      x$draws[, c(".draw", "time_h", "temp", "survival")]
    ))
  }
  stop("Could not recognize input as a predict_survival_curves() or ",
       "predict_heat_injury() result.", call. = FALSE)
}

# Internal: predict_survival_curves draws_matrix [ndraws x ngrid] → long.
surv_grid_to_long <- function(psc) {
  mat  <- psc$draws_matrix
  grid <- psc$grid
  if (ncol(mat) != nrow(grid))
    stop("draws_matrix and grid dimensions disagree.", call. = FALSE)
  ndraws <- nrow(mat)
  tibble::as_tibble(do.call(rbind, lapply(seq_len(ncol(mat)), function(j) {
    data.frame(.draw    = seq_len(ndraws),
               temp     = grid$temp[j],
               duration = grid$duration[j],
               survival = mat[, j])
  })))
}

# Internal: typeguard for extract_tdt() results.
stop_if_not_extract_tdt <- function(x) {
  if (!is.list(x) || is.null(x$z) || is.null(x$CTmax))
    stop("Expected an extract_tdt() result.", call. = FALSE)
  invisible(NULL)
}
