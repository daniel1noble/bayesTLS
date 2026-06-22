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
#' The `treedepth_max_attempted` field is the model's `max_treedepth` setting
#' (passed to `brms::brm()`); a saturation is a post-warmup transition where
#' the sampler hit that ceiling. Per-chain BFMI is computed from the energy
#' diagnostic in `brms::nuts_params(.)` following the standard Stan
#' definition (Var(ΔE)/Var(E)).
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
  td_vals <- subset(np, Parameter == "treedepth__")$Value
  td_max  <- if (length(td_vals)) max(td_vals, na.rm = TRUE) else NA_integer_
  treed   <- if (is.finite(td_max)) sum(td_vals >= td_max) else 0L

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
#' (so the result is parameterisation- and coding-agnostic \u2014 no coefficient-name
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
