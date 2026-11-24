# fit_lognormal(): thin exported wrapper around the internal engine tdt_fit()
# (R/fit_4pl.R) with the model fixed. Same formals as fit_4pl(); the curve
# itself is defined in R/models.R.

#' Fit the joint Bayesian four-parameter log-normal (probit) dose-response model
#'
#' Identical to [fit_4pl()] in every argument and in everything derived from
#' the fit, except that the curve is the log-normal (probit) one:
#' `p = low + (up - low) * Phi(-k * (log10(t) - mid))`, where `Phi` is the
#' standard normal distribution function. Individual death times are
#' log-normally distributed.
#'
#' `mid` is still the log10 time at which `p = (low + up) / 2`, so the TDT
#' line, the relative threshold, `z` and `CTmax` mean the same thing as in
#' [fit_4pl()]. The prior mean of `logk` is shifted so that the prior on the
#' slope of the curve at the midpoint matches the log-logistic's. See the
#' supplement's re-parameterisation table for the curve, its
#' absolute-threshold correction and its hazard.
#'
#' @inheritParams fit_4pl
#' @inherit fit_4pl return
#' @examples
#' \dontrun{
#' wf <- fit_lognormal(std)
#' tls(wf)
#' }
#' @export
fit_lognormal <- function(data,
                          random_effects = NULL,
                          temp_effects   = c("low", "up", "k", "mid"),
                          lower          = 0,
                          upper          = 1,
                          family         = NULL,
                          prior          = NULL,
                          chains         = 4,
                          iter           = 4000,
                          warmup         = floor(iter / 2),
                          cores          = chains,
                          seed           = 123,
                          backend        = "cmdstanr",
                          control        = list(adapt_delta = 0.95,
                                                max_treedepth = 12),
                          init           = 0,
                          file           = NULL,
                          file_refit     = "on_change",
                          fit            = TRUE,
                          bounds         = NULL,
                          ctmax          = NULL,
                          z              = NULL,
                          up             = NULL,
                          low            = NULL,
                          k              = NULL,
                          mid            = NULL,
                          by             = NULL,
                          threshold      = c("relative", "absolute"),
                          p              = 0.5,
                          t_ref          = 60,
                          data2          = NULL,
                          fix_asymptotes = NULL,
                          ...) {
  tdt_fit(model = "log-normal",
          data           = data,
          random_effects = random_effects,
          temp_effects   = temp_effects,
          lower          = lower,
          upper          = upper,
          family         = family,
          prior          = prior,
          chains         = chains,
          iter           = iter,
          warmup         = warmup,
          cores          = cores,
          seed           = seed,
          backend        = backend,
          control        = control,
          init           = init,
          file           = file,
          file_refit     = file_refit,
          fit            = fit,
          bounds         = bounds,
          ctmax          = ctmax,
          z              = z,
          up             = up,
          low            = low,
          k              = k,
          mid            = mid,
          by             = by,
          threshold      = threshold,
          p              = p,
          t_ref          = t_ref,
          data2          = data2,
          fix_asymptotes = fix_asymptotes,
          ...)
}
