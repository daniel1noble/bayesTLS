# Lightweight diagnostic helpers for a fitted tdt_4pl_workflow. Wraps the
# usual brms / posterior diagnostics into a single one-line summary plus a
# tidy parameter table, and points to brms::pp_check() for visual PPC.

#' Sampling diagnostics for a fitted TDT workflow
#'
#' Returns a tibble with the per-fit summary numbers a reviewer will want at
#' a glance: max Rhat, minimum bulk and tail ESS, divergent transitions, and
#' tree-depth saturations. Healthy values: Rhat < 1.01, ESS > 400, zero
#' divergences and treedepth hits.
#'
#' @param workflow A fitted `tdt_4pl_workflow`.
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

  divs  <- sum(subset(np, Parameter == "divergent__")$Value)
  treed <- sum(subset(np, Parameter == "treedepth__")$Value >=
               max(np$Value[np$Parameter == "treedepth__"], na.rm = TRUE))

  rh           <- brms::rhat(fit)
  rhat_max     <- max(rh, na.rm = TRUE)
  ess_summary  <- posterior::summarise_draws(
    posterior::as_draws(fit), "ess_bulk", "ess_tail"
  )
  ess_bulk_min <- min(ess_summary$ess_bulk, na.rm = TRUE)
  ess_tail_min <- min(ess_summary$ess_tail, na.rm = TRUE)

  tibble::tibble(
    rhat_max        = round(rhat_max, 4),
    ess_bulk_min    = round(ess_bulk_min),
    ess_tail_min    = round(ess_tail_min),
    divergences     = divs,
    treedepth_hits  = treed,
    rhat_pass       = rhat_max < 1.01,
    ess_pass        = ess_bulk_min > 400 & ess_tail_min > 400,
    divergence_pass = divs == 0
  )
}

#' Posterior parameter table on the natural scale
#'
#' Pulls the population-level posterior of the four reparameterised 4PL
#' parameters (`low`, `up`, `k`, `mid` intercept) and the `mid` temperature
#' slope, transformed back to the natural scale via the model's bounds. Also
#' includes `z = -1 / mid_temp`. Returns a one-row-per-parameter tibble with
#' median + 95% CrI.
#'
#' @param workflow A fitted `tdt_4pl_workflow`.
#' @return A tibble with columns `parameter`, `median`, `lower`, `upper`.
#' @examples
#' \dontrun{
#' tdt_parameter_table(wf)
#' }
#' @export
tdt_parameter_table <- function(workflow) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  d <- posterior::as_draws_df(workflow$fit) |> as.data.frame()
  b <- workflow$meta$bounds

  draws <- tibble::tibble(
    low      = b$low_min + stats::plogis(d$b_lowraw_Intercept) * b$low_w,
    up       = b$up_min  + stats::plogis(d$b_upraw_Intercept)  * b$up_w,
    k        = exp(d$b_logk_Intercept),
    mid_int  = d$b_mid_Intercept,
    mid_temp = d$b_mid_temp_c,
    z        = -1 / d$b_mid_temp_c
  )

  summary_row <- function(x, name) {
    q <- stats::quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE, names = FALSE)
    tibble::tibble(parameter = name,
                   median = q[2], lower = q[1], upper = q[3])
  }

  dplyr::bind_rows(
    summary_row(draws$low,      "low (lower asymptote)"),
    summary_row(draws$up,       "up (upper asymptote)"),
    summary_row(draws$k,        "k (slope)"),
    summary_row(draws$mid_int,  "mid intercept (at T_bar)"),
    summary_row(draws$mid_temp, "mid temp_c slope"),
    summary_row(draws$z,        "z (°C)")
  )
}
