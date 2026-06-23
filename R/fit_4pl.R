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
#'                       Ignored in the direct CTmax/z parameterisation.
#' @param bounds         Length-2 numeric `c(lower, upper)` giving the
#'                       response-scale asymptote interval. Defaults to
#'                       `c(lower, upper)`. Preferred over the separate
#'                       `lower`/`upper` (which it supersedes — see those).
#' @param ctmax,z        One-sided formulas for the **direct** parameterisation.
#'                       Supplying either switches `make_4pl_formula()` to fit
#'                       `CTmax` (as `CTmaxdev = CTmax - mean(temp)`) and `z`
#'                       (as `logz`) directly as nonlinear parameters, with `mid`
#'                       reconstructed from them. Carry the fixed and random
#'                       effects for tolerance and sensitivity, e.g.
#'                       `~ life_stage + (1 | Date)`. Omit `z` for a single
#'                       pooled `z` (`logz ~ 1`).
#' @param up,low,k       Optional one-sided formulas for the asymptote and
#'                       steepness sub-parameters. In **direct** mode, when
#'                       omitted they inherit the `ctmax`/`z` **fixed** effects
#'                       crossed with `temp_c` (random effects are not inherited);
#'                       intercept-only `ctmax`/`z` gives `~ temp_c`. In
#'                       **midpoint** mode they default to the `temp_effects`
#'                       structure (or `by`, see below); an explicit formula
#'                       always overrides.
#' @param mid            **Midpoint** mode only. Optional one-sided formula for
#'                       the midpoint sub-parameter (whose `temp_c` slope is
#'                       \eqn{-1/z}). Must carry `temp_c`. Defaults to the
#'                       `temp_effects`/`by` structure. Errors in direct mode
#'                       (where `mid` is reconstructed from `ctmax`/`z`).
#' @param by             **Midpoint** mode only. Optional character vector of
#'                       moderator column(s). A shortcut that fits `~ temp_c * by`
#'                       on **all four** sub-parameters (`low`/`up`/`k`/`mid`) at
#'                       once -- i.e. every aspect of the curve, plus the
#'                       midpoint's temperature slope (\eqn{-1/z}), varies by
#'                       group. `by` is the more specific instruction, so it
#'                       overrides `temp_effects`. Equivalent to writing the four
#'                       `~ temp_c * by` formulas by hand; use explicit formulas
#'                       for finer per-parameter control. Errors in direct mode
#'                       (group there via `ctmax = ~ 0 + <moderator>`).
#' @param threshold      `"relative"` (default) fits the midpoint backbone so
#'                       `CTmax`/`z` are the relative-threshold quantities;
#'                       `"absolute"` bakes the asymmetry correction
#'                       \eqn{C(T) = (1/k)\log((up-p)/(p-low))} into `mid` so
#'                       they are the absolute (`p`-survival) quantities.
#' @param p              Survival fraction defining the absolute threshold.
#'                       Default `0.5`.
#' @param log10_tref     `log10` of the reference exposure (in the model's time
#'                       unit) at which `CTmax` is defined. Default `0` (one time
#'                       unit). [fit_4pl()] derives this from `t_ref` and the
#'                       data's duration unit.
#' @return A `brmsformula` object.
#' @examples
#' make_4pl_formula()
#' make_4pl_formula(family = binomial(link = "identity"))
#' make_4pl_formula(temp_effects = "mid")            # constant-shape TDT
#' make_4pl_formula(family = brms::Beta(link = "identity"),
#'                  response_var = "survival")
#' make_4pl_formula(ctmax = ~ life_stage, z = ~ life_stage)  # direct CTmax/z
#' @export
make_4pl_formula <- function(random_effects = NULL,
                             bounds         = NULL,
                             family         = brms::beta_binomial(
                               link = "identity"),
                             response_var   = "survival",
                             temp_effects   = c("low", "up", "k", "mid"),
                             ctmax          = NULL,
                             z              = NULL,
                             up             = NULL,
                             low            = NULL,
                             k              = NULL,
                             mid            = NULL,
                             by             = NULL,
                             threshold      = c("relative", "absolute"),
                             p              = 0.5,
                             log10_tref     = 0,
                             lower          = 0,
                             upper          = 1) {

  if (is.null(bounds)) bounds <- c(lower, upper)
  threshold <- match.arg(threshold)
  b <- compute_4pl_bounds(bounds[1], bounds[2])

  low_expr <- sprintf("(%.6f + inv_logit(lowraw) * %.6f)", b$low_min, b$low_w)
  up_expr  <- sprintf("(%.6f + inv_logit(upraw)  * %.6f)", b$up_min,  b$up_w)

  # Response LHS + optional logit() wrap, shared across parameterisations.
  is_beta  <- identical(family$family, "beta")
  resp_lhs <- if (is_beta) response_var else "n_surv | trials(n_total)"
  wrap     <- function(rhs) if (is_beta && identical(family$link, "logit"))
                              sprintf("logit(%s)", rhs) else rhs

  if (is.null(ctmax) && is.null(z)) {
    ## ---- midpoint parameterisation ----
    temp_effects <- match.arg(temp_effects, c("low", "up", "k", "mid"),
                              several.ok = TRUE)
    main_rhs  <- sprintf(
      "%s + (%s - %s) / (1 + exp(exp(logk) * (logd - mid)))",
      low_expr, up_expr, low_expr)
    # Each sub-parameter's RHS. An explicit formula always wins. Otherwise:
    #  * with a `by` moderator -> ALL four get `temp_c * by` (treatment-coded, so
    #    the reference level keeps the centred Intercept prior and the per-group
    #    offsets are weakly shrunk toward it). `by` is the more specific
    #    instruction, so it overrides `temp_effects`.
    #  * without `by` -> the classic `temp_effects` default: `temp_c` for a
    #    sub-parameter in `temp_effects`, else `1` (constant in temperature).
    # Every sub-parameter carries at least `temp_c` in the all-four default; `mid`
    # MUST carry temp_c (its slope is -1/z, the core TDT quantity) -- enforced
    # below so e.g. temp_effects = c("low", "up") (mid dropped) errors clearly.
    by_term <- if (!is.null(by)) paste(by, collapse = " * ") else NULL
    par_rhs <- function(par) {
      if (!is.null(by_term)) return(paste0("temp_c * ", by_term))
      if (par %in% temp_effects) "temp_c" else "1"
    }
    res <- function(explicit, par) {
      if (!is.null(explicit)) formula_rhs(explicit, par_rhs(par)) else par_rhs(par)
    }
    low_r <- res(low, "low")
    up_r  <- res(up,  "up")
    k_r   <- res(k,   "k")
    mid_r <- res(mid, "mid")
    if (!grepl("temp_c", mid_r))
      stop("`mid` must always carry temp_c (its slope is -1/z, the core TDT ",
           "quantity); supply a `mid`/`by` formula that includes temp_c, or keep ",
           "\"mid\" in `temp_effects`.", call. = FALSE)
    # Guard the silent-pooling footgun: if low/up/k carry a moderator that mid
    # does NOT, z is pooled across that moderator even though the asymptotes vary
    # by it. That is a valid model, but rarely what a user writing
    # `low = ~ temp_c * species` intends -- point them at `by=` / `mid=`.
    shape_mods <- setdiff(rhs_fixed_vars(paste(low_r, up_r, k_r, sep = " + ")),
                          rhs_fixed_vars(mid_r))
    if (length(shape_mods))
      warning(sprintf(
        paste0("low/up/k vary by %s but mid does not, so z (= -1/mid's temp_c ",
               "slope) is POOLED across %s. Use by = \"%s\" (applies the moderator ",
               "to all four sub-parameters) or mid = ~ temp_c * %s to let z vary ",
               "by group."),
        paste(shape_mods, collapse = ", "), paste(shape_mods, collapse = ", "),
        shape_mods[1], shape_mods[1]), call. = FALSE)
    mid_rhs <- paste(c(mid_r, tdt_format_random_effects(random_effects)),
                     collapse = " + ")
    return(brms::bf(
      stats::as.formula(sprintf("%s ~ %s", resp_lhs, wrap(main_rhs))),
      stats::as.formula(paste("lowraw ~", low_r)),
      stats::as.formula(paste("upraw  ~", up_r)),
      stats::as.formula(paste("logk   ~", k_r)),
      stats::as.formula(paste("mid    ~", mid_rhs)),
      nl = TRUE, family = family))
  }

  ## ---- direct CTmax/z parameterisation ----
  # CTmaxdev and logz are estimated; low/up/k and mid are nlf-derived. mid is the
  # relative backbone, minus C(T) under threshold = "absolute".
  ctmax_rhs <- formula_rhs(ctmax, "1")
  z_rhs     <- formula_rhs(z,     "1")
  inherit   <- if (!is.null(ctmax)) ctmax_rhs else z_rhs
  # Absolute threshold builds log((up - p)/(p - low)); the disjoint asymptote
  # reparameterisation keeps low <= low_max and up >= up_min, so this is finite
  # for every draw iff low_max < p < up_min. With sub-unit response bounds (e.g.
  # PSII in [0.85, 1]) the achievable range never reaches p = 0.5, which would
  # give log() of a negative number (an opaque Stan failure). Catch it early with
  # an actionable message; relative thresholds are well defined for any bounds.
  if (identical(threshold, "absolute") && !(b$low_max < p && p < b$up_min))
    stop(sprintf(
      paste0("threshold = \"absolute\" with target p = %g is not achievable for ",
             "response bounds [%g, %g]: the disjoint-bounds reparameterisation ",
             "splits the asymptotes at the midpoint %.3f, so baking an absolute ",
             "%g%% correction into the fit is only safe for p near %.3f (here the ",
             "achievable range is roughly (%.3f, %.3f)). Either fit with ",
             "threshold = \"relative\" (the per-draw midpoint, defined for any ",
             "bounds), or -- to get an LT%g (or any LTx) -- fit with the relative ",
             "threshold and derive it afterwards via ",
             "extract_tdt()/tls(target_surv = %g), which inverts the fitted ",
             "surface post hoc and works for any p."),
      p, bounds[1], bounds[2], b$midpoint, 100 * p, b$midpoint,
      b$low_max, b$up_min, 100 * p, p),
      call. = FALSE)
  mid_rel   <- sprintf("(%g - (temp_c - CTmaxdev) / exp(logz))", log10_tref)
  mid_expr  <- if (identical(threshold, "relative")) mid_rel else
    sprintf("(%s - (1/exp(logk)) * log((up - %g)/(%g - low)))", mid_rel, p, p)
  main_rhs  <- "low + (up - low) / (1 + exp(exp(logk) * (logd - mid)))"

  brms::bf(
    stats::as.formula(sprintf("%s ~ %s", resp_lhs, wrap(main_rhs))),
    brms::nlf(stats::as.formula(paste("low ~", low_expr))),
    brms::nlf(stats::as.formula(paste("up  ~", up_expr))),
    brms::nlf(stats::as.formula(paste("mid ~", mid_expr))),
    stats::as.formula(paste("lowraw   ~", resolve_shape(low, inherit))),
    stats::as.formula(paste("upraw    ~", resolve_shape(up,  inherit))),
    stats::as.formula(paste("logk     ~", resolve_shape(k,   inherit))),
    stats::as.formula(paste("CTmaxdev ~", ctmax_rhs)),
    stats::as.formula(paste("logz     ~", z_rhs)),
    nl = TRUE, family = family)
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
#' @param bounds         Length-2 `c(lower, upper)` asymptote interval;
#'                       supersedes `lower`/`upper` when supplied.
#' @param ctmax,z,up,low,k One-sided formulas selecting the **direct** CTmax/z
#'                       parameterisation (see [make_4pl_formula()]). Supplying
#'                       `ctmax` and/or `z` fits CTmax and z directly; `up`/`low`/`k`
#'                       are inherited from them when omitted. `up`/`low`/`k` also
#'                       take explicit formulas in the **midpoint** parameterisation.
#' @param mid,by         **Midpoint** parameterisation only (see
#'                       [make_4pl_formula()]). `mid` is an explicit one-sided
#'                       formula for the midpoint (must carry `temp_c`); `by` is a
#'                       moderator column-name shortcut applying the moderator to
#'                       all four sub-parameters (`low`/`up`/`k`/`mid`) at once, so
#'                       a per-group CTmax/z fit needs only `by = "<moderator>"`.
#'                       Either records the moderator in `meta$group_vars`, so
#'                       [extract_tdt()] / [tls()] auto-derive per-group quantities
#'                       with `by = `. Both error in direct mode.
#' @param threshold      `"relative"` (default) or `"absolute"`; the threshold the
#'                       fitted CTmax/z refer to (direct mode only).
#' @param p              Survival fraction for the `"absolute"` threshold (direct
#'                       mode). Default `0.5` (LT50). Must lie within the response
#'                       bounds' achievable range, else `fit_4pl()` errors (use a
#'                       relative threshold for sub-unit-bounded responses).
#' @param t_ref          Reference exposure for CTmax, in **minutes** (default 60,
#'                       i.e. one hour). Converted to the model's log10-time scale
#'                       via the data's `duration_unit`.
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
                    ...) {

  threshold <- match.arg(threshold)
  if (!is.null(bounds)) { lower <- bounds[1]; upper <- bounds[2] }
  direct <- !is.null(ctmax) || !is.null(z)
  if (direct && !is.null(random_effects))
    stop("In the direct CTmax/z parameterisation, random effects go INSIDE the ",
         "`ctmax`/`z` formulas, e.g. ctmax = ~ 0 + grp + (1 | batch). The ",
         "`random_effects` argument is for the midpoint parameterisation and ",
         "would otherwise be silently ignored here.", call. = FALSE)
  if (direct && !is.null(mid))
    stop("`mid` is a midpoint-parameterisation argument; in the direct CTmax/z ",
         "parameterisation the midpoint is derived from ctmax/z and cannot be ",
         "set directly.", call. = FALSE)
  if (direct && !is.null(by))
    stop("`by` builds the midpoint per-group structure; in the direct CTmax/z ",
         "parameterisation, group per moderator via ctmax = ~ 0 + <moderator> ",
         "(z / up / low / k then inherit it).", call. = FALSE)

  # Make group coding irrelevant: a single-factor treatment term (`~ G` or
  # `~ 1 + G`) spans the same model as cell-means (`~ 0 + G`), but only
  # cell-means receives per-level CTmax/z priors -- under treatment coding the
  # mean-zero contrast prior inflates non-reference groups' marginal prior
  # variance (~sqrt(2)x), so the two codings would NOT give the same posterior.
  # Normalising single-factor treatment terms to cell-means makes the
  # parameterisations genuinely equivalent with correct priors. Multi-term,
  # interaction, and continuous moderators are left untouched.
  if (direct) {
    ctmax <- normalize_cellmeans(ctmax, data)
    z     <- normalize_cellmeans(z,     data)
    up    <- normalize_cellmeans(up,    data)
    low   <- normalize_cellmeans(low,   data)
    k     <- normalize_cellmeans(k,     data)
  }

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

  # Reference exposure for CTmax (direct mode): t_ref is in minutes; place it on
  # the model's log10-time scale. log10_tref = log10(t_ref / minutes-per-unit).
  # Use the SAME alias-aware unit mapping as extract_tdt()/derive_*(); a bespoke
  # switch() here silently fell back to tm = 1 for any non-canonical label
  # ("h", "hr", "Hours", ...), desynchronising log10_tref from the extractors.
  # The NULL-unit sentinel "model_units" keeps tm = 1 (no real time unit given).
  unit       <- meta$duration_unit %||% "model_units"
  tm         <- if (identical(unit, "model_units")) 1 else tdt_unit_to_minutes(unit)
  log10_tref <- log10(t_ref / tm)

  formula <- make_4pl_formula(random_effects = random_effects,
                              lower          = lower,
                              upper          = upper,
                              family         = family,
                              response_var   = response_var,
                              temp_effects   = temp_effects,
                              ctmax          = ctmax,
                              z              = z,
                              up             = up,
                              low            = low,
                              k              = k,
                              mid            = mid,
                              by             = by,
                              threshold      = threshold,
                              p              = p,
                              log10_tref     = log10_tref)

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
                             prior_phi      = prior_phi,
                             ctmax          = ctmax,
                             z              = z,
                             up             = up,
                             low            = low,
                             k              = k)
  }

  # Fixed-effect moderators (e.g. species, life_stage) the fit varies over:
  # direct mode reads them off ctmax/z/up/low/k; midpoint mode off up/low/k/mid
  # plus the `by` shortcut. `grouped` is coding-independent, so the group-aware
  # readers (extract_tdt(by=), tls(by=), ...) auto-derive per-group quantities
  # under BOTH `~ 0 + G` and `~ 1 + G` instead of silently returning the
  # reference level.
  group_vars <- if (direct) direct_group_vars(ctmax, z, up, low, k)
                else direct_group_vars(NULL, NULL, up, low, k, mid, by)

  meta_full <- utils::modifyList(meta, list(
    random_effects   = random_effects,
    temp_effects     = match.arg(temp_effects, c("low", "up", "k", "mid"),
                                 several.ok = TRUE),
    lower            = lower,
    upper            = upper,
    bounds           = compute_4pl_bounds(lower, upper),
    response_type    = response_type,
    family           = family$family,
    link             = family$link,
    parameterization = if (direct) "direct" else "midpoint",
    threshold        = threshold,
    t_ref            = t_ref,
    log10_tref       = log10_tref,
    group_vars       = group_vars,
    grouped          = length(group_vars) > 0
  ))

  workflow <- structure(
    list(fit = NULL, data = data, formula = formula,
         prior = prior, meta = meta_full),
    class = "bayes_tls"
  )

  if (!fit) return(workflow)

  # The default priors set a general `class = "b"` catch-all on CTmaxdev/logz
  # (needed only when CTmax/z carry a continuous moderator or interaction)
  # alongside the per-level / Intercept priors. For the common cell-means or
  # intercept-only fit every coefficient already has a specific prior, so brms
  # emits a harmless "global prior ... will not be used in the model" note for
  # the redundant catch-all. Muffle only that note (the fit is unaffected and
  # the user cannot act on it); all other warnings pass through untouched.
  workflow$fit <- withCallingHandlers(
    brms::brm(
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
    ),
    warning = function(w) {
      if (grepl("will not be used in the model", conditionMessage(w),
                fixed = TRUE))
        invokeRestart("muffleWarning")
    }
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

# --- internal: formula-resolution helpers for the direct CTmax/z mode -------
# (Unexported. Kept after the documented functions so the roxygen @export blocks
# above attach to make_4pl_formula()/fit_4pl(), not to these helpers.)

# Fixed-effect moderator variable names in a formula RHS string, excluding
# `temp_c` and random-effect grouping factors. Empty character vector => no
# moderator. Coding-independent: both `0 + species` (cell-means) and `species`
# (treatment) yield "species".
rhs_fixed_vars <- function(rhs) {
  if (is.null(rhs) || !nzchar(rhs)) return(character(0))
  rhs <- gsub("\\([^|]*\\|[^)]*\\)", "", rhs)          # drop (x | g) random terms
  rhs <- trimws(gsub("\\+\\s*$|^\\s*\\+", "", rhs))
  if (rhs %in% c("", "1", "0")) return(character(0))
  setdiff(all.vars(stats::as.formula(paste("~", rhs))), "temp_c")
}

# Fixed-effect moderator variable names across the parameterisation's formulas,
# plus any explicit `by` moderator. Direct mode reads ctmax/z/up/low/k; midpoint
# mode reads up/low/k/mid (+ by). Empty character vector => single-group fit.
direct_group_vars <- function(ctmax, z, up, low, k, mid = NULL, by = NULL) {
  fixed_vars <- function(f) {
    if (is.null(f) || !inherits(f, "formula")) return(character(0))
    rhs_fixed_vars(formula_rhs(f, "1"))
  }
  unique(c(by, unlist(lapply(list(ctmax, z, up, low, k, mid), fixed_vars))))
}

# Rewrite a single-factor TREATMENT-coded one-sided formula (`~ G` or `~ 1 + G`,
# optionally carrying random-effect terms) to CELL-MEANS (`~ 0 + G [+ (re)]`).
# Treatment and cell-means coding span the same model, but only cell-means draws
# per-level CTmax/z priors; this rewrite makes the two codings fit the identical
# model with correct priors. Returns `f` unchanged unless the FIXED part is
# exactly one factor/character variable with an implicit intercept -- continuous
# moderators, interactions, `I()` terms, multiple fixed terms, and formulas that
# already drop the intercept (`0 +` / `- 1`) are left as-is.
normalize_cellmeans <- function(f, data) {
  if (is.null(f) || !inherits(f, "formula")) return(f)
  rhs   <- formula_rhs(f, "1")
  re    <- regmatches(rhs, gregexpr("\\([^)]*\\|[^)]*\\)", rhs))[[1]]  # (x | g) terms
  fe    <- gsub("\\([^)]*\\|[^)]*\\)", "", rhs)                        # fixed part only
  terms <- trimws(strsplit(fe, "\\+")[[1]])
  terms <- terms[nzchar(terms)]
  if (any(grepl("(^|[^[:alnum:]_.])0([^[:alnum:]_.]|$)|-\\s*1", terms)))
    return(f)                                       # already no-intercept
  fixed <- setdiff(terms, "1")                      # drop an explicit intercept
  if (length(fixed) != 1L) return(f)                # one fixed term only
  v <- fixed
  if (!grepl("^[A-Za-z.][A-Za-z0-9._]*$", v)) return(f)   # bare variable only
  if (is.null(data) || !v %in% names(data)) return(f)
  if (!(is.factor(data[[v]]) || is.character(data[[v]]))) return(f)  # factors only
  stats::as.formula(
    paste("~", paste(c(paste0("0 + ", v), re), collapse = " + ")))
}

# RHS of a one-sided (or two-sided) formula as a string; `default` when NULL.
formula_rhs <- function(f, default) {
  if (is.null(f)) return(default)
  if (!inherits(f, "formula"))
    stop("formula arguments (ctmax, z, up, low, k) must be one-sided formulas ",
         "or NULL, e.g. `~ life_stage + (1 | Date)`.", call. = FALSE)
  paste(deparse(f[[length(f)]]), collapse = " ")
}

# Resolve a shape sub-parameter (up/low/k) RHS. An explicit formula wins.
# Otherwise inherit the CTmax/z FIXED effects crossed with temp_c (so the
# absolute correction can differ per group); RANDOM effects are NOT inherited
# (they stay on CTmax/z); intercept-only CTmax/z gives `temp_c`.
resolve_shape <- function(shape_f, inherit_rhs) {
  if (!is.null(shape_f)) return(formula_rhs(shape_f, "temp_c"))
  fe <- trimws(gsub("\\+?\\s*\\([^|]*\\|[^)]*\\)", "", inherit_rhs))  # drop (x | g)
  fe <- trimws(gsub("\\+\\s*$|^\\s*\\+\\s*", "", fe))
  if (fe %in% c("", "1")) return("temp_c")
  if (grepl("^0\\s*\\+", fe)) {                                       # cell-means 0 + G
    core <- trimws(sub("^0\\s*\\+\\s*", "", fe))
    if (grepl("\\+", core)) paste0("0 + ", core, " + temp_c:(", core, ")")
    else                    paste0("0 + ", core, " + temp_c:", core)
  } else {
    paste0("(", fe, ") * temp_c")                                    # treatment G * temp_c
  }
}
