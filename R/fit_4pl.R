# Build and fit the joint Bayesian 4PL with built-in beta_binomial(link =
# "identity"). The asymptote reparam keeps the 4PL output inside (low_min,
# up_max) so identity link is safe; see
# notes/2026-05-12-identity-vs-logit-link.qmd for the empirical comparison.

#' Build the brms formula for the joint Bayesian 4PL
#'
#' Constructs the joint 4PL with the disjoint-bounds reparameterisation on the
#' asymptotes (`low`, `up`) and a positivity reparam on `k`. All four 4PL
#' sub-parameters (`lowraw`, `upraw`, `logk`, `mid`) get a `temp_c` slope — the
#' model assumes temperature can affect every aspect of the dose-response curve.
#' Random intercepts (if any) attach to `mid`.
#'
#' The response term depends on the family:
#'
#' - **Count families** (`binomial`, `beta_binomial`): the response is
#'   `n_surv | trials(n_total) ~ <4PL>`.
#' - **`Beta` family** (continuous proportion): the response is
#'   `<response_var> ~ <4PL>` with no `trials()` term. With `link = "identity"`
#'   the 4PL value is the mean directly (safe because the reparam keeps it in
#'   `(0, 1)`); with `link = "logit"` the 4PL is wrapped in `logit()` so the mean
#'   is still the 4PL value.
#'
#' If a parameter has no real temperature effect, its `temp_c` slope shrinks
#' toward zero under the prior and the fit reduces cleanly to the simpler case.
#'
#' @param random_effects Optional character vector of grouping variables for
#'                       random intercepts on `mid`.
#' @param lower,upper    Response-scale bounds for the asymptotes. Default
#'                       `0, 1` for proportion data; set narrower for sublethal
#'                       data whose response sits well above `0`.
#' @param family         brms family. Default `beta_binomial(link = "identity")`.
#'                       Pass `binomial(link = "identity")` for the simpler
#'                       no-overdispersion case, or `brms::Beta(link = "identity")`
#'                       for a continuous proportion response.
#' @param response_var   Response column name for a `Beta` (continuous
#'                       proportion) fit. Ignored for count families. Default
#'                       `"survival"` (the column [standardize_data()] writes for
#'                       a `proportion` response).
#' @param temp_effects   Character vector naming which 4PL sub-parameters carry a
#'                       `temp_c` slope: any subset of `c("low", "up", "k",
#'                       "mid")`. Default is all four — temperature can affect
#'                       every aspect of the dose-response curve. `"mid"` must
#'                       always be present (its slope is \eqn{-1/z}, the core TDT
#'                       quantity). Restricting to `"mid"` gives the classical
#'                       *constant-shape* TDT model (asymptotes and slope shared
#'                       across temperature), which is the right choice for sparse
#'                       designs where the richer all-four model over-fits.
#' @return A `brmsformula` object.
#' @examples
#' make_4pl_formula()
#' make_4pl_formula(family = binomial(link = "identity"))
#' make_4pl_formula(temp_effects = "mid")            # constant-shape TDT
#' make_4pl_formula(family = brms::Beta(link = "identity"),
#'                  response_var = "survival")
#' @export
make_4pl_formula <- function(random_effects = NULL,
                             lower          = 0,
                             upper          = 1,
                             family         = brms::beta_binomial(
                               link = "identity"),
                             response_var   = "survival",
                             temp_effects   = c("low", "up", "k", "mid")) {

  temp_effects <- match.arg(temp_effects, c("low", "up", "k", "mid"),
                            several.ok = TRUE)
  if (!"mid" %in% temp_effects)
    stop("`mid` must always carry temp_c (its slope is -1/z, the core TDT ",
         "quantity); include \"mid\" in `temp_effects`.", call. = FALSE)

  b <- compute_4pl_bounds(lower, upper)

  low_expr <- sprintf("(%.6f + inv_logit(lowraw) * %.6f)", b$low_min, b$low_w)
  up_expr  <- sprintf("(%.6f + inv_logit(upraw)  * %.6f)", b$up_min,  b$up_w)

  main_rhs <- sprintf(
    "%s + (%s - %s) / (1 + exp(exp(logk) * (logd - mid)))",
    low_expr, up_expr, low_expr
  )

  # Each sub-parameter is "~ temp_c" if selected, else "~ 1" (constant in
  # temperature). Random intercepts always attach to mid.
  par_rhs   <- function(par) if (par %in% temp_effects) "temp_c" else "1"
  mid_terms <- c(par_rhs("mid"), tdt_format_random_effects(random_effects))
  mid_rhs   <- paste(mid_terms, collapse = " + ")

  family_name <- family$family
  if (identical(family_name, "beta")) {
    # Continuous proportion: no trials() term. Wrap in logit() only when the
    # family link is logit, so the mean is the 4PL value under either link.
    resp_rhs <- if (identical(family$link, "logit"))
                  sprintf("logit(%s)", main_rhs) else main_rhs
    response_formula <- sprintf("%s ~ %s", response_var, resp_rhs)
  } else {
    response_formula <- sprintf("n_surv | trials(n_total) ~ %s", main_rhs)
  }

  brms::bf(
    stats::as.formula(response_formula),
    stats::as.formula(paste("lowraw ~", par_rhs("low"))),
    stats::as.formula(paste("upraw  ~", par_rhs("up"))),
    stats::as.formula(paste("logk   ~", par_rhs("k"))),
    stats::as.formula(paste("mid    ~", mid_rhs)),
    nl = TRUE,
    family = family
  )
}

#' Fit the joint Bayesian 4PL to standardised TDT data
#'
#' Wraps [make_4pl_formula()] + [make_4pl_priors()] + [brms::brm()]. Returns
#' a workflow object containing the fit, data, formula, prior, and metadata
#' that downstream helpers (e.g. [extract_tdt()]) read from.
#'
#' By default all four 4PL sub-parameters get a `temp_c` slope, and random
#' intercepts attach to `mid` only. Use `temp_effects = "mid"` for the classical
#' constant-shape TDT model (recommended for sparse designs — see `temp_effects`).
#' To fit separate categories (life stages, species, populations), filter `data`
#' per category and call this function once per subset.
#'
#' @param data           Output of [standardize_data()].
#' @param random_effects Optional character vector of grouping variables for
#'                       random intercepts on `mid`.
#' @param temp_effects   Which 4PL sub-parameters carry a `temp_c` slope; a subset
#'                       of `c("low", "up", "k", "mid")` (default all four).
#'                       `"mid"` is always required. The all-four default needs
#'                       roughly >= 15 cells per fixed effect to be stable; for
#'                       sparse designs (few temperatures x durations, small n)
#'                       prefer `temp_effects = "mid"` (constant-shape), which
#'                       estimates only the midpoint's temperature slope (= -1/z)
#'                       and shares the asymptotes and steepness across
#'                       temperature. Passed to [make_4pl_formula()].
#' @param lower,upper    Response-scale bounds for the asymptotes. Default
#'                       `0, 1`; pass `lower = 0.85, upper = 1` for sublethal
#'                       data bounded above zero.
#' @param family         brms family. `NULL` (default) picks the family from the
#'                       data's `response_type` metadata:
#'                       `beta_binomial(link = "identity")` for count data and
#'                       `brms::Beta(link = "identity")` for a continuous
#'                       `proportion` response. Pass an explicit family to
#'                       override (e.g. `binomial(link = "identity")` for the
#'                       no-overdispersion count case).
#' @param prior          Optional `brmsprior` object. If `NULL` (default),
#'                       [make_4pl_priors()] builds defaults from `data`,
#'                       `lower`, `upper`, and `random_effects`.
#' @param chains,iter,warmup,cores,seed Sampling arguments passed to [brms::brm()].
#' @param backend        Either `"cmdstanr"` (default) or `"rstan"`.
#' @param control        List of HMC control parameters.
#' @param init           Initial values for the sampler. Default `0` initialises
#'                       every unconstrained parameter at zero.
#' @param file           Optional path (without extension) to cache the fit.
#' @param file_refit     Passed to [brms::brm()]. Default `"on_change"`.
#' @param fit            Logical. If `FALSE`, returns the workflow spec without
#'                       fitting — useful for inspecting the formula and priors.
#' @param ...            Further arguments passed to [brms::brm()].
#' @return A list of class `"bayes_tls"` with elements `fit`, `data`,
#'         `formula`, `prior`, `meta`. See [print.bayes_tls()],
#'         [summary.bayes_tls()], and [plot.bayes_tls()] for the available
#'         methods.
#' @examples
#' \dontrun{
#' d  <- standardize_data(raw, ...)
#' wf <- fit_4pl(d, fit = FALSE)         # inspect spec only
#' wf <- fit_4pl(d, file = "output/models/my_fit")
#' wf <- fit_4pl(d, lower = 0.85)        # sublethal data
#' }
#' @export
fit_4pl <- function(data,
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
                    ...) {

  meta <- attr(data, "tdt_meta") %||%
    list(temp_mean      = mean(data$temp, na.rm = TRUE),
         duration_unit  = "model_units",
         random_effects = random_effects,
         response_type  = "count")

  response_type <- meta$response_type %||% "count"
  response_var  <- meta$response_var %||% "survival"

  # Default family from the response type when not supplied.
  if (is.null(family)) {
    family <- if (identical(response_type, "proportion"))
                brms::Beta(link = "identity")
              else
                brms::beta_binomial(link = "identity")
  }

  formula <- make_4pl_formula(random_effects = random_effects,
                              lower          = lower,
                              upper          = upper,
                              family         = family,
                              response_var   = response_var,
                              temp_effects   = temp_effects)

  if (is.null(prior)) {
    # Set the phi prior only for families that carry a precision parameter
    # (beta and beta_binomial); skip it for binomial (no overdispersion).
    family_name <- family$family
    prior_phi   <- if (family_name %in% c("beta_binomial", "beta"))
                     "gamma(2, 0.1)" else NULL
    prior <- make_4pl_priors(data           = data,
                             lower          = lower,
                             upper          = upper,
                             random_effects = random_effects,
                             prior_phi      = prior_phi)
  }

  meta_full <- utils::modifyList(meta, list(
    random_effects = random_effects,
    temp_effects   = match.arg(temp_effects, c("low", "up", "k", "mid"),
                               several.ok = TRUE),
    lower          = lower,
    upper          = upper,
    bounds         = compute_4pl_bounds(lower, upper),
    response_type  = response_type,
    family         = family$family,
    link           = family$link
  ))

  workflow <- structure(
    list(fit = NULL, data = data, formula = formula,
         prior = prior, meta = meta_full),
    class = "bayes_tls"
  )

  if (!fit) return(workflow)

  workflow$fit <- brms::brm(
    formula    = formula,
    data       = data,
    prior      = prior,
    chains     = chains,
    iter       = iter,
    warmup     = warmup,
    cores      = cores,
    seed       = seed,
    backend    = backend,
    control    = control,
    init       = init,
    file       = file,
    file_refit = file_refit,
    ...
  )

  workflow
}

#' Whether a workflow has been fitted
#'
#' @param workflow Object returned by [fit_4pl()].
#' @return Logical scalar.
#' @examples
#' wf <- list(fit = NULL); class(wf) <- "bayes_tls"
#' has_fit(wf)
#' @export
has_fit <- function(workflow) {
  inherits(workflow, "bayes_tls") && !is.null(workflow$fit)
}

#' Extract the fitted brms model from a workflow
#'
#' Returns the underlying [brms::brmsfit] object held in a `bayes_tls`
#' workflow, ready for any brms / posterior / bayesplot helper (e.g.
#' [brms::bayes_R2()], [brms::mcmc_plot()], [posterior::as_draws_df()]).
#' Prefer this over reaching into `workflow$fit` directly: it errors clearly
#' on an unfitted workflow instead of silently returning `NULL`. Pairs with
#' the predicate [has_fit()].
#'
#' @param workflow A fitted `bayes_tls` workflow returned by [fit_4pl()].
#' @return The `brmsfit` object stored in the workflow.
#' @examples
#' \dontrun{
#' wf  <- fit_4pl(std)
#' fit <- get_brmsfit(wf)
#' brms::bayes_R2(fit)
#' }
#' @export
get_brmsfit <- function(workflow) {
  if (!has_fit(workflow))
    stop("workflow has no fit; call fit_4pl() first.", call. = FALSE)
  workflow$fit
}
