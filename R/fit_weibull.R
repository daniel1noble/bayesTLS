# fit_weibull(): thin exported wrapper around the internal engine tdt_fit()
# (R/fit_4pl.R) with the model fixed. Same formals as fit_4pl(); the curve
# itself is defined in R/models.R.

#' Fit the joint Bayesian four-parameter Weibull dose-response model
#'
#' Identical to [fit_4pl()] in every argument and in everything derived from
#' the fit, except for the curve.
#'
#' * `type = 1` (default), Weibull type 1:
#'   `p = low + (up - low) * 2^(-exp(k * (log10(t) - mid)))`.
#'   This is the Gompertz-shaped curve on log time; individual death times
#'   are Weibull-distributed.
#' * `type = 2`, Weibull type 2:
#'   `p = up - (up - low) * 2^(-exp(-k * (log10(t) - mid)))`.
#'   Individual death times are Frechet-distributed (inverse Weibull).
#'
#' Both are written with base 2 so that `mid` is the log10 time at which
#' `p = (low + up) / 2`, exactly as in [fit_4pl()]; the TDT line, the relative
#' threshold, `z` and `CTmax` therefore mean the same thing. The standard
#' (inflection-point) location of the curve sits `0.367 / k` above (type 1) or
#' below (type 2) `mid`. The prior mean of `logk` is shifted so that the prior
#' on the slope at the midpoint matches the log-logistic's.
#'
#' @param type `1` (default) for Weibull type 1, `2` for Weibull type 2.
#' @inheritParams fit_4pl
#' @inherit fit_4pl return
#' @examples
#' \dontrun{
#' wf <- fit_weibull(std)              # type 1
#' wf2 <- fit_weibull(std, type = 2)
#' }
#' @export
fit_weibull <- function(data,
                        type           = 1,
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
  if (!(is.numeric(type) && length(type) == 1L && type %in% c(1, 2)))
    stop("`type` must be 1 (Weibull type 1, Gompertz-shaped) or 2 (Weibull ",
         "type 2).", call. = FALSE)
  tdt_fit(model = if (type == 1) "weibull1" else "weibull2",
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
