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
#' Returns a list with the model spec, a posterior summary of every sampled
#' parameter (population-level coefficients `b_*`, random-effect SDs `sd_*`,
#' and overdispersion `phi` if present), and headline HMC diagnostics (max
#' Rhat, min ESS, divergent transitions, treedepth saturations).
#'
#' For natural-scale 4PL parameters (`low`, `up`, `k`, `z`), use
#' [tdt_parameter_table()]. For the TDT quantities (`z`, `CTmax_1hr`,
#' optionally `T_crit`), use [extract_tdt()].
#'
#' @param object A fitted `bayes_tls` workflow.
#' @param ...    Ignored.
#' @return A list with class `"summary.bayes_tls"` containing `meta`,
#'         `shape`, `parameters` (a `posterior::draws_summary` tibble), and
#'         `diagnostics` (named list of scalar HMC summaries).
#' @examples
#' \dontrun{
#' s <- summary(wf)
#' s$diagnostics
#' s$parameters
#' }
#' @export
summary.bayes_tls <- function(object, ...) {
  if (!has_fit(object))
    stop("workflow has no fit; call fit_4pl() first.", call. = FALSE)

  fit   <- object$fit
  draws <- posterior::as_draws(fit)
  params <- posterior::summarise_draws(draws)

  np      <- brms::nuts_params(fit)
  divs    <- sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
  td_vals <- np$Value[np$Parameter == "treedepth__"]
  td_max  <- if (length(td_vals)) max(td_vals, na.rm = TRUE) else NA_real_
  treed   <- if (is.finite(td_max)) sum(td_vals >= td_max) else 0L

  structure(
    list(
      meta = object$meta,
      shape = list(
        n_obs       = nrow(object$data),
        n_temps     = length(unique(object$data$temp)),
        n_durations = length(unique(object$data$duration)),
        n_draws     = brms::ndraws(fit)
      ),
      parameters  = params,
      diagnostics = list(
        max_rhat       = if ("rhat"     %in% names(params))
                           max(params$rhat,     na.rm = TRUE) else NA_real_,
        min_ess_bulk   = if ("ess_bulk" %in% names(params))
                           min(params$ess_bulk, na.rm = TRUE) else NA_real_,
        min_ess_tail   = if ("ess_tail" %in% names(params))
                           min(params$ess_tail, na.rm = TRUE) else NA_real_,
        divergences    = as.integer(divs),
        treedepth_hits = as.integer(treed)
      )
    ),
    class = "summary.bayes_tls"
  )
}

#' Print method for `summary.bayes_tls`
#'
#' @param x      A `summary.bayes_tls` object.
#' @param digits Digits for the posterior summary table. Default `3`.
#' @param ...    Ignored.
#' @return The object, invisibly.
#' @export
print.summary.bayes_tls <- function(x, digits = 3, ...) {
  bounds <- x$meta$bounds
  cat("<bayes_tls summary>\n")
  cat("  Data:    ", x$shape$n_obs,        "rows;",
      x$shape$n_temps,     "temperatures;",
      x$shape$n_durations, "durations\n")
  cat("  T_bar:   ", round(x$meta$temp_mean, 2), "\n")
  cat(sprintf(
    "  Bounds:   response in (%.3f, %.3f); low in (%.3f, %.3f); up in (%.3f, %.3f)\n",
    x$meta$lower, x$meta$upper,
    bounds$low_min, bounds$low_max,
    bounds$up_min,  bounds$up_max
  ))
  if (!is.null(x$meta$random_effects))
    cat("  RE:      ", paste(x$meta$random_effects, collapse = ", "), "\n")
  cat("  Draws:   ", x$shape$n_draws, "\n\n")

  d <- x$diagnostics
  cat("HMC diagnostics:\n")
  cat(sprintf("  max Rhat = %.4f   min ESS bulk = %s   min ESS tail = %s\n",
              d$max_rhat,
              if (is.finite(d$min_ess_bulk)) format(round(d$min_ess_bulk)) else "NA",
              if (is.finite(d$min_ess_tail)) format(round(d$min_ess_tail)) else "NA"))
  cat(sprintf("  divergences = %d   treedepth hits = %d\n\n",
              d$divergences, d$treedepth_hits))

  cat("Posterior summary:\n")
  print(as.data.frame(x$parameters), row.names = FALSE, digits = digits)
  invisible(x)
}

#' MCMC trace plot for a fitted `bayes_tls` workflow
#'
#' Post-warmup trace of each sampled parameter, chains coloured. Use to
#' eyeball chain mixing alongside the numeric Rhat / ESS in
#' [summary.bayes_tls()].
#'
#' @param x    A fitted `bayes_tls` workflow.
#' @param pars Optional character vector of parameter names to plot. Default
#'             `NULL` selects population-level coefficients (`b_*`),
#'             random-effect SDs (`sd_*`), and the dispersion parameter
#'             (`phi`) when the family carries one.
#' @param ...  Ignored.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' plot(wf)
#' plot(wf, pars = c("b_mid_Intercept", "b_mid_temp_c"))
#' }
#' @export
plot.bayes_tls <- function(x, pars = NULL, ...) {
  if (!has_fit(x))
    stop("workflow has no fit; call fit_4pl() first.", call. = FALSE)

  draws    <- posterior::as_draws_df(x$fit)
  all_pars <- setdiff(names(draws), c(".chain", ".iteration", ".draw"))

  if (is.null(pars)) {
    pars <- grep("^(b_|sd_|phi$)", all_pars, value = TRUE)
    if (length(pars) == 0L) pars <- all_pars
  } else {
    missing_pars <- setdiff(pars, all_pars)
    if (length(missing_pars))
      stop("Unknown parameter(s): ", paste(missing_pars, collapse = ", "),
           call. = FALSE)
  }

  long <- do.call(rbind, lapply(pars, function(p) {
    data.frame(parameter = p,
               iteration = draws$.iteration,
               chain     = factor(draws$.chain),
               value     = draws[[p]])
  }))

  ggplot2::ggplot(long,
                  ggplot2::aes(x = iteration, y = value, colour = chain)) +
    ggplot2::geom_line(alpha = 0.7, linewidth = 0.3) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y") +
    theme_tdt() +
    ggplot2::labs(x = "iteration (post-warmup)", y = NULL, colour = "Chain")
}
