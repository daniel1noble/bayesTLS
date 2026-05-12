# Build and fit the joint Bayesian 4PL with built-in beta_binomial(link =
# "identity"). The asymptote reparam keeps the 4PL output inside (low_min,
# up_max) so identity link is safe; see
# notes/2026-05-12-identity-vs-logit-link.qmd for the empirical comparison.

#' Build the brms formula for the joint Bayesian 4PL
#'
#' Constructs `n_surv | trials(n_total) ~ <4PL>` with the disjoint-bounds
#' reparameterisation on the asymptotes (`low`, `up`) and a positivity reparam
#' on `k`. All four 4PL sub-parameters (`lowraw`, `upraw`, `logk`, `mid`) get
#' a `temp_c` slope — the model assumes temperature can affect every aspect of
#' the dose-response curve. Random intercepts (if any) attach to `mid`. Family
#' is `beta_binomial(link = "identity")`.
#'
#' If a parameter has no real temperature effect, its `temp_c` slope shrinks
#' toward zero under the prior and the fit reduces cleanly to the simpler case.
#'
#' @param random_effects Optional character vector of grouping variables for
#'                       random intercepts on `mid`.
#' @param lower,upper    Response-scale bounds for the asymptotes. Default
#'                       `0, 1` for proportion data; set narrower for sublethal
#'                       data whose response sits well above `0`.
#' @return A `brmsformula` object.
#' @examples
#' make_4pl_formula()
#' make_4pl_formula(random_effects = c("Date", "Tank"))
#' make_4pl_formula(lower = 0.85, upper = 1)   # PSII-like
#' @export
make_4pl_formula <- function(random_effects = NULL,
                             lower          = 0,
                             upper          = 1) {

  b <- compute_4pl_bounds(lower, upper)

  low_expr <- sprintf("(%.6f + inv_logit(lowraw) * %.6f)", b$low_min, b$low_w)
  up_expr  <- sprintf("(%.6f + inv_logit(upraw)  * %.6f)", b$up_min,  b$up_w)

  main_rhs <- sprintf(
    "%s + (%s - %s) / (1 + exp(exp(logk) * (logd - mid)))",
    low_expr, up_expr, low_expr
  )

  mid_terms <- c("temp_c", tdt_format_random_effects(random_effects))
  mid_rhs   <- paste(mid_terms, collapse = " + ")

  brms::bf(
    stats::as.formula(sprintf("n_surv | trials(n_total) ~ %s", main_rhs)),
    stats::as.formula("lowraw ~ temp_c"),
    stats::as.formula("upraw  ~ temp_c"),
    stats::as.formula("logk   ~ temp_c"),
    stats::as.formula(paste("mid    ~", mid_rhs)),
    nl = TRUE,
    family = brms::beta_binomial(link = "identity")
  )
}

#' Fit the joint Bayesian 4PL to standardised TDT data
#'
#' Wraps [make_4pl_formula()] + [make_4pl_priors()] + [brms::brm()]. Returns
#' a workflow object containing the fit, data, formula, prior, and metadata
#' that downstream helpers (e.g. [extract_tdt()]) read from.
#'
#' The model is fixed: all four 4PL sub-parameters get a `temp_c` slope, and
#' random intercepts attach to `mid` only. To fit separate categories (life
#' stages, species, populations), filter `data` per category and call this
#' function once per subset.
#'
#' @param data           Output of [standardize_data()].
#' @param random_effects Optional character vector of grouping variables for
#'                       random intercepts on `mid`.
#' @param lower,upper    Response-scale bounds for the asymptotes. Default
#'                       `0, 1`; pass `lower = 0.85, upper = 1` for sublethal
#'                       data bounded above zero.
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
#' @return A list of class `"tdt_4pl_workflow"` with elements `fit`, `data`,
#'         `formula`, `prior`, `meta`.
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
                    lower          = 0,
                    upper          = 1,
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
         random_effects = random_effects)

  formula <- make_4pl_formula(random_effects = random_effects,
                              lower          = lower,
                              upper          = upper)

  if (is.null(prior)) {
    prior <- make_4pl_priors(data           = data,
                             lower          = lower,
                             upper          = upper,
                             random_effects = random_effects)
  }

  meta_full <- utils::modifyList(meta, list(
    random_effects = random_effects,
    lower          = lower,
    upper          = upper,
    bounds         = compute_4pl_bounds(lower, upper)
  ))

  workflow <- structure(
    list(fit = NULL, data = data, formula = formula,
         prior = prior, meta = meta_full),
    class = "tdt_4pl_workflow"
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
#' wf <- list(fit = NULL); class(wf) <- "tdt_4pl_workflow"
#' has_fit(wf)
#' @export
has_fit <- function(workflow) {
  inherits(workflow, "tdt_4pl_workflow") && !is.null(workflow$fit)
}

#' Brief printout of a TDT workflow
#'
#' @param x A `tdt_4pl_workflow` object.
#' @param ... Ignored.
#' @export
print.tdt_4pl_workflow <- function(x, ...) {
  bounds <- x$meta$bounds
  cat("<tdt_4pl_workflow>\n")
  cat("  Data:    ", nrow(x$data), "rows;",
      length(unique(x$data$temp)), "temperatures;",
      length(unique(x$data$duration)), "durations\n")
  cat("  T_bar:   ", round(x$meta$temp_mean, 2), "\n")
  cat(sprintf("  Range:    response in (%.3f, %.3f); low in (%.3f, %.3f); up in (%.3f, %.3f)\n",
              x$meta$lower, x$meta$upper,
              bounds$low_min, bounds$low_max,
              bounds$up_min,  bounds$up_max))
  if (!is.null(x$meta$random_effects))
    cat("  RE:      ", paste(x$meta$random_effects, collapse = ", "), "\n")
  cat("  Status:  ", if (is.null(x$fit)) "spec only (not fitted)"
                      else sprintf("fitted (%d draws)", brms::ndraws(x$fit)),
      "\n")
  invisible(x)
}
