#' Default weakly informative priors for the joint Bayesian 4PL
#'
#' Priors are weakly informative and do not adapt to the data (which would
#' understate uncertainty). They are constructed to match
#' [make_4pl_formula()]'s disjoint-bounds reparameterisation
#' `low = low_min + inv_logit(lowraw) * low_w` and
#' `up = up_min + inv_logit(upraw) * up_w` for the asymptote intervals derived
#' from `lower` and `upper` by [compute_4pl_bounds()].
#'
#' Each of the four 4PL sub-parameters (`lowraw`, `upraw`, `logk`, `mid`) gets
#' an Intercept-specific prior plus a general `class = "b"` prior covering the
#' `temp_c` slope. Random-effect SDs on `mid` get one `prior_random_sd` row per
#' grouping variable.
#'
#' @param data            Output of [standardize_data()]. Only `logd` is used,
#'                        to centre the midpoint Intercept prior.
#' @param lower,upper     Response-scale bounds for the asymptotes. Default
#'                        `0, 1` for proportion data. For sublethal data
#'                        bounded above `0` (e.g. photosystem-II between 0.85
#'                        and 1), set `lower` accordingly.
#' @param random_effects  Optional character vector of grouping variables for
#'                        random intercepts on `mid`. Adds one
#'                        `prior_random_sd` row per group.
#' @param prior_random_sd Stan prior string for random-effect SDs. Default
#'                        `"exponential(2)"`.
#' @param prior_phi       Stan prior string for the beta-binomial precision
#'                        parameter `phi`. Default `"gamma(2, 0.1)"`. Pass
#'                        `NULL` to omit the `phi` prior — needed when fitting
#'                        with `family = binomial()` (no overdispersion).
#' @return A `brmsprior` object ready to pass to [fit_4pl()].
#' @examples
#' raw <- data.frame(
#'   temperature_C = rep(c(30, 32, 34), each = 4),
#'   exposure_h    = rep(c(1, 2, 4, 8), times = 3),
#'   n             = 30L,
#'   alive         = c(29, 28, 25, 5, 30, 27, 18, 2, 28, 22, 10, 1)
#' )
#' d <- standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
#'                       n_total = "n", n_surv = "alive")
#' make_4pl_priors(d)
#' @export
make_4pl_priors <- function(data,
                            lower           = 0,
                            upper           = 1,
                            random_effects  = NULL,
                            prior_random_sd = "exponential(2)",
                            prior_phi       = "gamma(2, 0.1)") {

  b <- compute_4pl_bounds(lower, upper)

  mid_start <- stats::median(data$logd, na.rm = TRUE)

  # Centre the raw parameters so low ~ lower + 0.02 * range and
  # up ~ lower + 0.98 * range — the typical TDT asymptotes for any
  # response-scale range.
  low_target  <- lower + 0.02 * (upper - lower)
  up_target   <- lower + 0.98 * (upper - lower)
  lowraw_mean <- stats::qlogis((low_target - b$low_min) / b$low_w)
  upraw_mean  <- stats::qlogis((up_target  - b$up_min)  / b$up_w)
  logk_mean   <- log(2)

  priors <- c(
    brms::set_prior(sprintf("normal(%.6f, 1)", lowraw_mean),
                    class = "b", nlpar = "lowraw", coef = "Intercept"),
    brms::set_prior("normal(0, 0.5)",
                    class = "b", nlpar = "lowraw"),

    brms::set_prior(sprintf("normal(%.6f, 1)", upraw_mean),
                    class = "b", nlpar = "upraw", coef = "Intercept"),
    brms::set_prior("normal(0, 0.5)",
                    class = "b", nlpar = "upraw"),

    brms::set_prior(sprintf("normal(%.6f, 1)", logk_mean),
                    class = "b", nlpar = "logk", coef = "Intercept"),
    brms::set_prior("normal(0, 0.3)",
                    class = "b", nlpar = "logk"),

    brms::set_prior(sprintf("normal(%.6f, 1.5)", mid_start),
                    class = "b", nlpar = "mid", coef = "Intercept"),
    brms::set_prior("normal(0, 0.6)",
                    class = "b", nlpar = "mid")
  )

  if (!is.null(prior_phi)) {
    priors <- c(priors, brms::set_prior(prior_phi, class = "phi"))
  }

  for (re_var in tdt_random_effect_variables(random_effects)) {
    priors <- c(priors,
                brms::set_prior(prior_random_sd,
                                class = "sd", nlpar = "mid", group = re_var))
  }

  priors
}
