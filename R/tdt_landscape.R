# Two-dimensional survival surface across temperature × duration. Same
# machinery as predict_survival_curves(), but defaults to a denser 2D grid
# suitable for a heatmap / contour visualisation.

#' Posterior survival landscape across the temperature × duration plane
#'
#' Returns the posterior median (and 95% credible band) of survival probability
#' on a 2-D grid covering the experimental range. Suitable for plotting as a
#' heatmap or contour with [plot_tdt_landscape()] (Phase 1d). For curves at a
#' handful of temperatures, use [predict_survival_curves()] instead.
#'
#' @param workflow      Fitted `bayes_tls`.
#' @param temp_grid     Numeric vector of temperatures (°C). Default: 120
#'                      equally spaced values across the training-data range.
#' @param duration_grid Numeric vector of durations. Default: 120 equally
#'                      spaced values across the training-data range. A linear
#'                      (rather than log-spaced) default keeps the grid regular
#'                      on the linear duration axis that [plot_tdt_landscape()]
#'                      uses by default, so the survival surface renders without
#'                      uneven-raster artefacts.
#' @param ndraws        Posterior draws to use; `NULL` (default) uses the full
#'                      posterior. Pass an integer to subsample for speed.
#' @param probs         Quantile probabilities. Default `c(0.025, 0.5, 0.975)`.
#' @return A list with the same shape as [predict_survival_curves()] (a
#'         `summary` tibble with `temp`, `duration`, `survival_median`,
#'         `survival_lower`, `survival_upper`, plus `draws_matrix` and `grid`).
#' @examples
#' \dontrun{
#' wf  <- fit_4pl(d, ...)
#' lsp <- derive_tdt_landscape(wf)
#' lsp$summary
#' }
#' @export
derive_tdt_landscape <- function(workflow,
                                 temp_grid     = NULL,
                                 duration_grid = NULL,
                                 ndraws        = NULL,
                                 probs         = c(0.025, 0.5, 0.975)) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  data <- workflow$data
  if (is.null(temp_grid)) {
    trange    <- range(data$temp, na.rm = TRUE)
    temp_grid <- seq(trange[1], trange[2], length.out = 120)
  }
  if (is.null(duration_grid)) {
    drange        <- range(data$duration, na.rm = TRUE)
    duration_grid <- seq(drange[1], drange[2], length.out = 120)
  }

  predict_survival_curves(
    workflow  = workflow,
    temps     = temp_grid,
    durations = duration_grid,
    ndraws    = ndraws,
    probs     = probs
  )
}
