# Lightweight diagnostic helpers for a fitted bayes_tls. Wraps the
# usual brms / posterior diagnostics into a single one-line summary plus a
# tidy parameter table, and points to brms::pp_check() for visual PPC.

#' Sampling diagnostics for a fitted TDT workflow
#'
#' Returns a tibble with the per-fit summary numbers a reviewer will want at
#' a glance: max Rhat, minimum bulk and tail ESS, divergent transitions,
#' tree-depth saturations, minimum BFMI per chain, and pass flags for each
#' criterion plus a combined `all_pass`. Healthy values:
#' \itemize{
#'   \item Rhat < 1.01
#'   \item ESS bulk and tail > 400
#'   \item zero divergent transitions
#'   \item zero tree-depth saturations
#'   \item BFMI > 0.3 per chain
#' }
#'
#' The `treedepth_max` field is the model's `max_treedepth` setting (the ceiling
#' passed to `brms::brm()`, default 10); a saturation is a post-warmup transition
#' that hit that ceiling. A fit that merely tops out *below* the ceiling has zero
#' saturations. Per-chain BFMI is computed from the energy diagnostic in
#' `brms::nuts_params(.)` following the standard Stan definition (Var(ΔE)/Var(E)).
#'
#' @param workflow A fitted `bayes_tls`.
#' @return A tibble with one row of diagnostic statistics.
#' @examples
#' \dontrun{
#' diagnose_tdt_fit(wf)
#' }
#' @export
diagnose_tdt_fit <- function(workflow) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  fit <- workflow$fit
  np  <- brms::nuts_params(fit)

  divs    <- sum(subset(np, Parameter == "divergent__")$Value)
  td_vals    <- subset(np, Parameter == "treedepth__")$Value
  td_ceiling <- tdt_max_treedepth(fit)
  # A saturation hits the configured CEILING. An iteration merely sitting at the
  # run's observed maximum (when that maximum is below the ceiling) is NOT a
  # saturation -- the previous `sum(td_vals >= max(td_vals))` flagged healthy fits.
  treed      <- if (length(td_vals)) sum(td_vals >= td_ceiling, na.rm = TRUE) else 0L

  rh           <- brms::rhat(fit)
  rhat_max     <- max(rh, na.rm = TRUE)
  ess_summary  <- posterior::summarise_draws(
    posterior::as_draws(fit), "ess_bulk", "ess_tail"
  )
  ess_bulk_min <- min(ess_summary$ess_bulk, na.rm = TRUE)
  ess_tail_min <- min(ess_summary$ess_tail, na.rm = TRUE)

  # BFMI per chain. nuts_params returns one row per (Chain, Iteration,
  # Parameter); "energy__" is the per-iteration energy. BFMI ≈ Var(ΔE) /
  # Var(E) per chain (post-warmup only — nuts_params is already post-warmup).
  energy <- subset(np, Parameter == "energy__")
  bfmi_per_chain <- if (nrow(energy) > 0L) {
    vapply(split(energy$Value, energy$Chain), function(e) {
      if (length(e) < 2L) return(NA_real_)
      dE <- diff(e)
      var_e <- stats::var(e)
      if (!is.finite(var_e) || var_e <= 0) NA_real_ else stats::var(dE) / var_e
    }, numeric(1))
  } else NA_real_
  bfmi_min <- if (all(is.na(bfmi_per_chain))) NA_real_
              else min(bfmi_per_chain, na.rm = TRUE)

  rhat_pass       <- isTRUE(rhat_max < 1.01)
  ess_pass        <- isTRUE(ess_bulk_min > 400 && ess_tail_min > 400)
  divergence_pass <- isTRUE(divs == 0L)
  treedepth_pass  <- isTRUE(treed == 0L)
  bfmi_pass       <- isTRUE(is.finite(bfmi_min) && bfmi_min > 0.3)
  all_pass <- rhat_pass && ess_pass && divergence_pass &&
              treedepth_pass && bfmi_pass

  tibble::tibble(
    rhat_max        = round(rhat_max, 4),
    ess_bulk_min    = round(ess_bulk_min),
    ess_tail_min    = round(ess_tail_min),
    divergences     = as.integer(divs),
    treedepth_max   = as.integer(td_ceiling),
    treedepth_hits  = as.integer(treed),
    bfmi_min        = round(bfmi_min, 4),
    rhat_pass       = rhat_pass,
    ess_pass        = ess_pass,
    divergence_pass = divergence_pass,
    treedepth_pass  = treedepth_pass,
    bfmi_pass       = bfmi_pass,
    all_pass        = all_pass
  )
}

# Configured max_treedepth ceiling for a fitted brms model, across backends.
# cmdstanr exposes it in the run metadata; rstan stores it in the per-chain
# stan_args control list; if neither is reachable we fall back to Stan's default
# of 10. Returns a single numeric (the ceiling, max across chains if they differ).
tdt_max_treedepth <- function(fit) {
  bf <- fit$fit
  md <- tryCatch(bf$metadata(), error = function(e) NULL)            # cmdstanr
  if (!is.null(md) && !is.null(md$max_treedepth)) {
    mt <- suppressWarnings(as.numeric(md$max_treedepth))
    mt <- mt[is.finite(mt)]
    if (length(mt)) return(max(mt))
  }
  sa <- tryCatch(methods::slot(bf, "stan_args"), error = function(e) NULL)  # rstan
  if (!is.null(sa)) {
    mt <- suppressWarnings(as.numeric(unlist(lapply(sa, function(a) a$control$max_treedepth))))
    mt <- mt[is.finite(mt)]
    if (length(mt)) return(max(mt))
  }
  10
}

#' Posterior parameter table on the natural scale
#'
#' Pulls the population-level posterior of the reparameterised 4PL parameters,
#' transformed back to the natural scale via the model's bounds, as a
#' one-row-per-parameter tibble with median + 95% CrI.
#'
#' For the **midpoint** parameterisation the rows are `low`, `up`, `k`, the
#' `mid` intercept and temperature slope, and `z = -1 / mid_temp`. For the
#' **direct** CTmax/z parameterisation the rows are `low`, `up`, `k`,
#' `CTmax` (at the model's reference dose) and `z`. All quantities are read by
#' evaluating `brms::posterior_linpred(nlpar=)` at `temp_c = 0` and `temp_c = 1`
#' (so the result is parameterisation- and coding-agnostic -- no coefficient-name
#' parsing). For a fit whose CTmax/z vary by a moderator the table is returned
#' per group, with the moderator column(s) prepended.
#'
#' @param workflow A fitted `bayes_tls`.
#' @param by Optional character vector of moderator columns to report per group.
#'   `NULL` (default) uses the fit's moderators (`meta$group_vars`); a
#'   single-condition fit then returns one block (no group column).
#' @return A tibble with columns `parameter`, `median`, `lower`, `upper`
#'   (plus the moderator column(s) for a grouped fit).
#' @examples
#' \dontrun{
#' tdt_parameter_table(wf)
#' }
#' @export
tdt_parameter_table <- function(workflow, by = NULL) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  fit    <- get_brmsfit(workflow)
  by     <- tdt_resolve_by(workflow, by)
  direct <- identical(workflow$meta$parameterization, "direct")
  l10    <- workflow$meta$log10_tref %||% 0
  tbar   <- workflow$meta$temp_mean  %||% 0

  # Evaluate low/up/k/mid at temp_c = 0 (constant-shape asymptotes + midpoint
  # intercept) and temp_c = 1 (for the linear midpoint slope), per group.
  nd <- tls_build_grid(fit$data, by = by, temp = "temp_c", temp_grid = c(0, 1))
  sp <- tls_eval_subpars(fit, nd, workflow$meta$bounds, mode = "relative")

  summary_row <- function(x, name) {
    # Non-finite per-draw values (e.g. +/-Inf z/CTmax from a near-zero midpoint
    # slope) are set to NA before quantiling, matching tls_local_z() in
    # R/tls_engine.R. Otherwise na.rm = TRUE would propagate +/-Inf into the
    # reported bounds and could corrupt the median.
    x[!is.finite(x)] <- NA_real_
    q <- stats::quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE, names = FALSE)
    tibble::tibble(parameter = name, median = q[2], lower = q[1], upper = q[3])
  }

  per_group <- lapply(unique(nd$.grp), function(g) {
    ci <- which(nd$.grp == g)
    c0 <- ci[which(nd$temp_c[ci] == 0)]; c1 <- ci[which(nd$temp_c[ci] == 1)]
    low <- sp$low[, c0]; up <- sp$up[, c0]; k <- sp$k[, c0]
    mid_int <- sp$mid[, c0]; mid_temp <- sp$mid[, c1] - sp$mid[, c0]
    z <- -1 / mid_temp
    tab <- if (direct) dplyr::bind_rows(
      summary_row(low, "low (lower asymptote)"),
      summary_row(up,  "up (upper asymptote)"),
      summary_row(k,   "k (slope)"),
      summary_row(tbar + (l10 - mid_int) / mid_temp, "CTmax (\u00b0C, at reference dose)"),
      summary_row(z,   "z (\u00b0C)")
    ) else dplyr::bind_rows(
      summary_row(low,      "low (lower asymptote)"),
      summary_row(up,       "up (upper asymptote)"),
      summary_row(k,        "k (slope)"),
      summary_row(mid_int,  "mid intercept (at T_bar)"),
      summary_row(mid_temp, "mid temp_c slope"),
      summary_row(z,        "z (\u00b0C)")
    )
    if (!is.null(by))
      tab <- cbind(nd[ci[1], by, drop = FALSE], tab, row.names = NULL)
    tab
  })
  dplyr::bind_rows(per_group)
}

#' Bayesian \eqn{R^2} for a fitted TDT workflow
#'
#' A thin tidy wrapper around [brms::bayes_R2()] for a fitted `bayes_tls`
#' workflow: it pulls the underlying `brmsfit` (erroring clearly on an unfitted
#' workflow, via [get_brmsfit()]) and returns the posterior median \eqn{R^2}
#' with its standard error and 95% credible interval as a one-row tibble, rather
#' than the bare matrix `brms::bayes_R2()` returns. Pass `...` straight through
#' to [brms::bayes_R2()] (e.g. `re_formula = NA` to exclude group-level effects,
#' or `ndraws =` to subsample). Map over a named list of workflows to build a
#' multi-fit table.
#'
#' @param workflow A fitted `bayes_tls` workflow returned by [fit_4pl()].
#' @param ... Further arguments passed to [brms::bayes_R2()].
#' @return A one-row tibble with columns `estimate`, `est_error`, `lower`, and
#'   `upper` (the lower and upper credible bounds, 2.5% and 97.5% by default;
#'   controlled by `probs` passed via `...`).
#' @seealso [brms::bayes_R2()], [diagnose_tdt_fit()]
#' @examples
#' \dontrun{
#' wf <- fit_4pl(std)
#' bayes_R2_tls(wf)
#' # Multiple fits in one table:
#' purrr::imap_dfr(list(binom = wf1, beta = wf2),
#'                 ~ cbind(model = .y, bayes_R2_tls(.x)))
#' }
#' @export
bayes_R2_tls <- function(workflow, ...) {
  fit <- get_brmsfit(workflow)            # clear error if the workflow is unfitted
  r   <- brms::bayes_R2(fit, ...)
  # brms::bayes_R2() returns a 1-row matrix with columns
  # Estimate / Est.Error / Q2.5 / Q97.5 (the quantile labels track `probs`).
  # This wrapper reshapes that summary; passing summary = FALSE via ... yields a
  # raw draws matrix without those columns, which this function cannot summarise.
  if (is.null(dim(r)) || !"Estimate" %in% colnames(r))
    stop("bayes_R2_tls() needs the summarised output of brms::bayes_R2(); ",
         "do not pass summary = FALSE.", call. = FALSE)
  cn  <- colnames(r)
  qlo <- grep("^Q", cn)[1]
  qhi <- grep("^Q", cn)[length(grep("^Q", cn))]
  tibble::tibble(
    estimate  = unname(r[1, "Estimate"]),
    est_error = unname(r[1, "Est.Error"]),
    lower     = unname(r[1, qlo]),
    upper     = unname(r[1, qhi])
  )
}
