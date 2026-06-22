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
#' @param bounds          Length-2 `c(lower, upper)` asymptote interval;
#'                        supersedes `lower`/`upper` when supplied.
#' @param ctmax,z,up,low,k One-sided formulas selecting the **direct** CTmax/z
#'                        priors (matching [make_4pl_formula()]). Supplying
#'                        `ctmax`/`z` switches to direct-mode priors: the same
#'                        centred asymptote priors (per factor level for
#'                        cell-means terms), weakly-informative `CTmaxdev`/`logz`
#'                        priors, and random-effect SD priors on `CTmaxdev`/`logz`
#'                        only (never on the shape sub-parameters). `up`/`low`/`k`
#'                        follow the same inheritance resolution as the formula.
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
                            prior_phi       = "gamma(2, 0.1)",
                            bounds          = NULL,
                            ctmax           = NULL,
                            z               = NULL,
                            up              = NULL,
                            low             = NULL,
                            k               = NULL) {

  if (is.null(bounds)) bounds <- c(lower, upper)
  b <- compute_4pl_bounds(bounds[1], bounds[2])

  # Centre the raw asymptotes so low ~ bounds[1] + 0.02 * range and
  # up ~ bounds[1] + 0.98 * range — the typical TDT asymptotes for any range.
  low_target  <- bounds[1] + 0.02 * (bounds[2] - bounds[1])
  up_target   <- bounds[1] + 0.98 * (bounds[2] - bounds[1])
  lowraw_mean <- stats::qlogis((low_target - b$low_min) / b$low_w)
  upraw_mean  <- stats::qlogis((up_target  - b$up_min)  / b$up_w)
  logk_mean   <- log(2)

  if (is.null(ctmax) && is.null(z)) {
    ## ---- midpoint priors (unchanged) ----
    mid_start <- stats::median(data$logd, na.rm = TRUE)
    priors <- c(
      brms::set_prior(sprintf("normal(%.6f, 1)", lowraw_mean),
                      class = "b", nlpar = "lowraw", coef = "Intercept"),
      brms::set_prior("normal(0, 0.5)", class = "b", nlpar = "lowraw"),
      brms::set_prior(sprintf("normal(%.6f, 1)", upraw_mean),
                      class = "b", nlpar = "upraw", coef = "Intercept"),
      brms::set_prior("normal(0, 0.5)", class = "b", nlpar = "upraw"),
      brms::set_prior(sprintf("normal(%.6f, 1)", logk_mean),
                      class = "b", nlpar = "logk", coef = "Intercept"),
      brms::set_prior("normal(0, 0.3)", class = "b", nlpar = "logk"),
      brms::set_prior(sprintf("normal(%.6f, 1.5)", mid_start),
                      class = "b", nlpar = "mid", coef = "Intercept"),
      brms::set_prior("normal(0, 0.6)", class = "b", nlpar = "mid")
    )
    if (!is.null(prior_phi))
      priors <- c(priors, brms::set_prior(prior_phi, class = "phi"))
    for (re_var in tdt_random_effect_variables(random_effects))
      priors <- c(priors, brms::set_prior(prior_random_sd,
                                          class = "sd", nlpar = "mid", group = re_var))
    return(priors)
  }

  ## ---- direct CTmax/z priors ----
  # Centred asymptote priors are load-bearing (a flat prior gave 191 divergences;
  # see notes/2026-06-21-prior-sensitivity-direct-ctmax-z.qmd). They match the
  # midpoint centres exactly. CTmax/z are well identified and take simple priors;
  # random effects attach to CTmaxdev/logz only (never inherited onto the shape).
  ctmax_rhs <- formula_rhs(ctmax, "1")
  z_rhs     <- formula_rhs(z,     "1")
  inherit   <- if (!is.null(ctmax)) ctmax_rhs else z_rhs
  shape_rhs <- list(lowraw = resolve_shape(low, inherit),
                    upraw  = resolve_shape(up,  inherit),
                    logk   = resolve_shape(k,   inherit))

  asy <- function(centre, nlpar, slope_sd) {
    cells <- cell_means_coefs(shape_rhs[[nlpar]], data)
    centred <- if (is.null(cells))
      list(brms::set_prior(sprintf("normal(%.6f, 1)", centre),
                           class = "b", nlpar = nlpar, coef = "Intercept"))
    else lapply(cells, function(cf)
      brms::set_prior(sprintf("normal(%.6f, 1)", centre),
                      class = "b", nlpar = nlpar, coef = cf))
    c(do.call(c, centred),
      brms::set_prior(sprintf("normal(0, %s)", slope_sd), class = "b", nlpar = nlpar))
  }

  priors <- c(
    asy(lowraw_mean, "lowraw", "0.5"),
    asy(upraw_mean,  "upraw",  "0.5"),
    asy(logk_mean,   "logk",   "0.3"),
    brms::set_prior("normal(0, 10)", class = "b", nlpar = "CTmaxdev"),
    brms::set_prior(sprintf("normal(%.6f, 0.7)", log(3)),
                    class = "b", nlpar = "logz")
  )
  if (!is.null(prior_phi))
    priors <- c(priors, brms::set_prior(prior_phi, class = "phi"))
  for (g in re_groups(ctmax))
    priors <- c(priors, brms::set_prior(prior_random_sd,
                                        class = "sd", nlpar = "CTmaxdev", group = g))
  for (g in re_groups(z))
    priors <- c(priors, brms::set_prior(prior_random_sd,
                                        class = "sd", nlpar = "logz", group = g))
  priors
}

# --- internal: structure helpers for the direct CTmax/z priors -------------
# (Unexported. Placed after make_4pl_priors() so its roxygen @export attaches to
# the function, not to these helpers.)

# Cell-means coefficient names for a "0 + G ..." resolved RHS (e.g.
# "0 + life_stage + temp_c:life_stage" -> c("life_stageA", "life_stageB", ...));
# NULL when the term is not cell-means coded (then the prior uses "Intercept").
cell_means_coefs <- function(rhs, data) {
  if (!grepl("^0\\s*\\+", rhs)) return(NULL)
  core <- trimws(sub("^0\\s*\\+\\s*", "", rhs))
  fac  <- trimws(strsplit(core, "[:+*]")[[1]][1])
  if (is.null(data) || !fac %in% names(data)) return(NULL)
  paste0(fac, levels(factor(data[[fac]])))
}

# Grouping-factor names in a one-sided formula's (x | g) random-effect terms.
re_groups <- function(f) {
  if (is.null(f)) return(character(0))
  rhs <- formula_rhs(f, "")
  m <- regmatches(rhs, gregexpr("\\|\\s*[^)]+", rhs))[[1]]
  trimws(sub("^\\|\\s*", "", m))
}
