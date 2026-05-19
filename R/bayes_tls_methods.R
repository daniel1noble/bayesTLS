# S3 methods for objects of class `bayes_tls` (currently produced by
# fit_4pl()). Provides print, summary, and plot (trace) methods so users can
# inspect a fit at three levels of detail: header, parameter+diagnostic
# table, MCMC mixing plot. The underlying brmsfit lives in `$fit` and brms /
# bayesplot methods are still available on it (e.g. `plot(wf$fit)`).

#' Compact print method for a `bayes_tls` workflow
#'
#' One-screen header: data shape, asymptote/response bounds, random-effect
#' grouping (if any), and fit status (draws if fitted, "spec only" otherwise).
#' Call [summary()][summary.bayes_tls()] for parameter posteriors and HMC
#' diagnostics, [plot()][plot.bayes_tls()] for MCMC trace plots.
#'
#' @param x   A `bayes_tls` object returned by [fit_4pl()].
#' @param ... Ignored.
#' @return The object, invisibly.
#' @examples
#' \dontrun{
#' wf <- fit_4pl(std)
#' print(wf)
#' }
#' @export
print.bayes_tls <- function(x, ...) {
  bounds <- x$meta$bounds
  cat("<bayes_tls>\n")
  cat("  Data:    ", nrow(x$data), "rows;",
      length(unique(x$data$temp)),     "temperatures;",
      length(unique(x$data$duration)), "durations\n")
  cat("  T_bar:   ", round(x$meta$temp_mean, 2), "\n")
  cat(sprintf(
    "  Bounds:   response in (%.3f, %.3f); low in (%.3f, %.3f); up in (%.3f, %.3f)\n",
    x$meta$lower, x$meta$upper,
    bounds$low_min, bounds$low_max,
    bounds$up_min,  bounds$up_max
  ))
  if (!is.null(x$meta$random_effects))
    cat("  RE:      ", paste(x$meta$random_effects, collapse = ", "), "\n")
  cat("  Status:  ",
      if (is.null(x$fit)) "spec only (not fitted)"
      else sprintf("fitted (%d draws)", brms::ndraws(x$fit)),
      "\n")
  invisible(x)
}

#' Summarise a fitted `bayes_tls` workflow
#'
#' Delegates to `summary()` on the underlying `brmsfit` (`x$fit`), so the
#' returned object is a `summary.brmsfit` with the population-level
#' coefficient table, group-level standard deviations, family parameters,
#' and HMC diagnostics already laid out by brms.
#'
#' The high-level workflow context (data shape, $\bar T$, asymptote
#' bounds, random-effect grouping, draw count) is available via
#' [print.bayes_tls()]. For natural-scale 4PL parameters (`low`, `up`,
#' `k`, `z`), use [tdt_parameter_table()]. For the TDT quantities
#' (`z`, `CTmax_1hr`, optionally `T_crit`), use [extract_tdt()].
#'
#' @param object A fitted `bayes_tls` workflow.
#' @param ...    Passed through to [brms::summary.brmsfit()] (e.g.
#'               `prob`, `mc_se`, `priors`, `robust`).
#' @return A `summary.brmsfit` object (brms handles printing).
#' @examples
#' \dontrun{
#' summary(wf)
#' summary(wf, prob = 0.9, robust = TRUE)
#' }
#' @export
summary.bayes_tls <- function(object, ...) {
  if (!has_fit(object))
    stop("workflow has no fit; call fit_4pl() first.", call. = FALSE)
  summary(object$fit, ...)
}

#' MCMC mixing plot for a fitted `bayes_tls` workflow
#'
#' Delegates to `plot()` on the underlying `brmsfit` (`x$fit`), which
#' produces brms's default `mcmc_combo` layout (per-parameter density
#' on the left, post-warmup trace on the right, chains coloured). Use
#' this to eyeball chain mixing alongside the numeric Rhat / ESS in
#' [summary.bayes_tls()].
#'
#' @param x    A fitted `bayes_tls` workflow.
#' @param ...  Passed through to `brms`'s `plot.brmsfit` method (e.g.
#'             `pars`, `combo`, `N`, `ask`).
#' @return The brms plot output (invisibly), typically a list of
#'         `bayesplot` ggplots.
#' @examples
#' \dontrun{
#' plot(wf)
#' plot(wf, pars = "^b_mid")
#' plot(wf, combo = c("dens_overlay", "trace"))
#' }
#' @export
plot.bayes_tls <- function(x, ...) {
  if (!has_fit(x))
    stop("workflow has no fit; call fit_4pl() first.", call. = FALSE)
  plot(x$fit, ...)
}
