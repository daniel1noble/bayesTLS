# Helper functions for the two-stage bias simulation.
#
# These functions are NOT part of the user-facing bayesTLS R package — they
# are simulation-study internals private to this manuscript. They depend on
# bayesTLS's public API (standardize_data, fit_4pl, extract_tdt,
# diagnose_tdt_fit) but the simulation harness itself is project-only.
#
# Source this file from the simulation driver:
#   source(here::here("scripts", "sim_twostage_helpers.R"))
#
# Roxygen tags (@param, @return, @export) are kept for readability only —
# they are not processed by devtools::document() since this file lives in
# scripts/ rather than R/.

# ----- truth (data-generating parameters; locked) -----------------------------
#
# Beta-binomial 4PL. Temperature acts on the midpoint always (slope m_beta1);
# DGP variants extend this by giving u, ell, and/or k linear temperature
# slopes via the `u_beta1`, `ell_beta1`, `k_beta1` fields of the truth list
# (defaulting to 0 for the baseline DGP). The analytical z and CTmax_1hr
# *targets* both estimators are judged against are the OLS slope/intercept of
# log10(LT50_{p=0.5})(T) at the simulation design temperatures — computed
# from the (T-varying) truth surface, not from a single-T linearisation. For
# the baseline DGP this OLS target reduces in closed form to the existing
# -1/m_beta1 and the existing CTmax_1hr formula.

#' Default truth parameters for the two-stage bias simulation
#'
#' Beta-binomial 4PL. The `dgp` argument selects a data-generating regime:
#'
#' \itemize{
#'   \item `"baseline"` — temperature acts only on the midpoint. u, ell, k constant in T.
#'   \item `"sym_ul"` — u and ell shift with T symmetrically; `(u+ell)/2` preserved.
#'   \item `"asym_u"` — u shifts with T, ell constant; midpoint displaced.
#'   \item `"varying_k"` — k shifts with T; u, ell constant.
#' }
#'
#' Returns a named list with the truth fields plus `dgp`, plus the OLS
#' z / CTmax_1hr targets the simulation scores against (see header comment in
#' this file for the truth-definition note).
#'
#' @param dgp Character. One of `"baseline"`, `"sym_ul"`, `"asym_u"`,
#'            `"varying_k"`.
#' @param u_0 Optional numeric — overrides the upper-asymptote intercept (the
#'            value of u at T = T_bar). `NULL` keeps the DGP preset.
#' @param u_beta1,ell_beta1,k_beta1 Optional numeric — override the temperature
#'            slopes of u, ell, k. `NULL` keeps the DGP preset.
#' @param design Character — `"full"` (5 temps × 6 durations, the original
#'            design) or `"sparse"` (3 temps × 4 durations). Controls the
#'            grid on which the OLS truth is evaluated.
#' @return Named list of truth parameters.
#' @export
sim_twostage_truth <- function(dgp       = "baseline",
                                u_0       = NULL,
                                u_beta1   = NULL,
                                ell_beta1 = NULL,
                                k_beta1   = NULL,
                                design    = "full") {
  base <- list(
    ell       = 0.05,
    u         = 0.92,
    k         = 8,
    m_beta0   = 1.5,
    m_beta1   = -0.15,
    T_bar     = 34,
    phi       = 5,
    u_beta1   = 0,    # slope of u   in (T - T_bar)
    ell_beta1 = 0,    # slope of ell in (T - T_bar)
    k_beta1   = 0     # slope of k   in (T - T_bar)
  )
  out <- switch(dgp,
    baseline  = base,
    sym_ul    = utils::modifyList(base, list(u_beta1 = -0.01,
                                              ell_beta1 = 0.01)),
    asym_u    = utils::modifyList(base, list(u_beta1 = -0.01)),
    varying_k = utils::modifyList(base, list(k_beta1 = 0.25)),
    stop("Unknown dgp: '", dgp, "'. ",
         "Use one of: baseline, sym_ul, asym_u, varying_k.")
  )

  # Explicit overrides take precedence over the DGP preset.
  if (!is.null(u_0))       out$u         <- u_0
  if (!is.null(u_beta1))   out$u_beta1   <- u_beta1
  if (!is.null(ell_beta1)) out$ell_beta1 <- ell_beta1
  if (!is.null(k_beta1))   out$k_beta1   <- k_beta1

  out$dgp    <- dgp
  out$design <- design

  # OLS targets — slope of log10(LT50_{p=0.5})(T) at design temperatures.
  tt <- compute_ols_truth(out, design = design)
  out$z_true         <- tt$z_true
  out$CTmax_1hr_true <- tt$CTmax_1hr_true
  out
}

#' Internal — OLS targets for z and CTmax_1hr at the design temperatures
#'
#' For each temperature in the design grid, compute
#' `log10(LT50_{p=0.5})(T) = mid_t(T) + (1/k(T)) * log((u(T)-0.5)/(0.5-ell(T)))`
#' using the (possibly T-varying) truth surface; then OLS-fit a line to those
#' values and read off `z_true = -1/slope` and `CTmax_1hr_true =
#' (log10(60) - intercept)/slope`. Used by [sim_twostage_truth()].
#'
#' @param p Named list of truth parameters (output of [sim_twostage_truth()]
#'          before `z_true`/`CTmax_1hr_true` are written).
#' @return List with `z_true` and `CTmax_1hr_true`.
#' @keywords internal
compute_ols_truth <- function(p, design = "full") {
  Ts <- switch(design,
               full   = c(30, 32, 34, 36, 38),
               sparse = c(30, 34, 38),
               stop("Unknown design: '", design,
                    "'. Use 'full' or 'sparse'.", call. = FALSE))
  T_c       <- Ts - p$T_bar
  u_T       <- p$u   + p$u_beta1   * T_c
  ell_T     <- p$ell + p$ell_beta1 * T_c
  k_T       <- p$k   + p$k_beta1   * T_c
  mid_T     <- p$m_beta0 + p$m_beta1 * T_c

  if (any(u_T <= 0.5) || any(ell_T >= 0.5) || any(k_T <= 0))
    stop("DGP gives degenerate u/ell/k at design temperatures: ",
         "need u > 0.5, ell < 0.5, k > 0.", call. = FALSE)

  ratio       <- (u_T - 0.5) / (0.5 - ell_T)
  log10_lt50  <- mid_T + (1 / k_T) * log(ratio)

  fit       <- stats::lm(log10_lt50 ~ Ts)
  co        <- stats::coef(fit)
  slope     <- co[["Ts"]]
  intercept <- co[["(Intercept)"]]

  list(
    z_true         = -1 / slope,
    CTmax_1hr_true = (log10(60) - intercept) / slope
  )
}

#' Default factorial design grid for the two-stage bias simulation
#'
#' Five assay temperatures and six durations matching the supplement's
#' simulation. Returns the (temp, duration) crossing only — replication is
#' added by [sim_twostage_dataset()] per scenario.
#'
#' @param design Character. `"full"` (5 temperatures × 6 durations) or
#'               `"sparse"` (3 temperatures × 4 durations). Default `"full"`.
#' @return Tibble with columns `T`, `t`, `log10_t`, `T_c`.
#' @export
sim_twostage_grid <- function(design = "full") {
  truth <- sim_twostage_truth()
  cells <- switch(design,
    full = list(T = c(30, 32, 34, 36, 38),
                t = c(1, 5, 15, 45, 135, 405)),
    sparse = list(T = c(30, 34, 38),
                  t = c(1, 15, 135, 405)),
    stop("Unknown design: '", design,
         "'. Use 'full' or 'sparse'.", call. = FALSE)
  )
  tidyr::expand_grid(T = cells$T, t = cells$t) |>
    dplyr::mutate(log10_t = log10(t),
                  T_c     = T - truth$T_bar)
}

#' Simulate one beta-binomial 4PL dataset
#'
#' Generates `n_reps` replicate cups per (temperature × duration) cell, each
#' with `N` individuals drawn uniformly from `n_ind_range`. Survival counts are
#' drawn from a beta-binomial with the truth's $p$ and dispersion $\phi$.
#'
#' @param n_reps      Replicate cups per cell.
#' @param n_ind_range Length-2 integer vector — discrete-uniform range for
#'                    `N` per cup. Default `c(10, 20)`.
#' @param seed        Integer seed for reproducibility.
#' @param truth       Optional truth-parameter list (defaults to
#'                    [sim_twostage_truth()]).
#' @return Tibble with one row per cup: `T`, `t`, `log10_t`, `T_c`, `rep`,
#'         `n`, `y` (alive count), `p_true`.
#' @export
sim_twostage_dataset <- function(n_reps,
                                 n_ind_range = c(10, 20),
                                 seed,
                                 truth = sim_twostage_truth()) {
  set.seed(seed)
  # The truth list carries the design label (set inside sim_twostage_truth);
  # use it so the (temp, duration) grid matches the design the OLS truth was
  # computed against.
  grid_design <- if (!is.null(truth$design)) truth$design else "full"
  grid   <- sim_twostage_grid(design = grid_design)
  design <- tidyr::expand_grid(grid, rep = seq_len(n_reps))

  # T-varying 4PL parameters (any beta1 = 0 reduces to the constant case).
  u_T    <- truth$u   + truth$u_beta1   * design$T_c
  ell_T  <- truth$ell + truth$ell_beta1 * design$T_c
  k_T    <- truth$k   + truth$k_beta1   * design$T_c
  mid_T  <- truth$m_beta0 + truth$m_beta1 * design$T_c

  p_true <- ell_T + (u_T - ell_T) /
    (1 + exp(k_T * (design$log10_t - mid_T)))

  n      <- sample(seq(n_ind_range[1], n_ind_range[2]),
                   size = nrow(design), replace = TRUE)
  alpha  <- p_true * truth$phi
  beta   <- (1 - p_true) * truth$phi
  p_draw <- stats::rbeta(nrow(design), shape1 = alpha, shape2 = beta)
  y      <- stats::rbinom(nrow(design), size = n, prob = p_draw)

  design$n      <- n
  design$y      <- y
  design$p_true <- p_true
  tibble::as_tibble(design)
}

# ----- classical two-stage pipeline -------------------------------------------

#' Classical two-stage TDT estimator (logit-GLM per temperature + OLS)
#'
#' Stage 1: for each temperature, fits a binomial GLM with logit link on
#' survival ~ log10(duration); reads off $\log_{10}\text{LT50}$ at $\hat p = 0.5$.
#' Stage 2: OLS of $\log_{10}\text{LT50}$ on temperature; reports $z = -1/\hat\beta_1$
#' and $CT_{max_{1hr}} = (\log_{10} 60 - \hat\beta_0)/\hat\beta_1$.
#'
#' 95% CIs come from the delta method on the Stage-2 coefficients. This is the
#' implementation field practice uses; alternative weighting (1/SE^2) at
#' Stage 2 is intentionally not used here so the simulation benchmarks the
#' pipeline as published.
#'
#' @param data      Output of [sim_twostage_dataset()] (raw simulated dataset).
#' @param t_ref_min Reference exposure time, in minutes. Default 60.
#' @return Named list. `$success` is `FALSE` if Stage 1 failed for any
#'         temperature (degenerate cells or non-finite LT50). Otherwise:
#'         `$z` and `$CTmax_1hr` each as `list(point, lower, upper)`;
#'         `$stage1` (per-temperature LT50 estimates).
#' @export
fit_two_stage_classical <- function(data, t_ref_min = 60) {
  temps <- sort(unique(data$T))

  per_temp <- lapply(temps, function(T_value) {
    d <- subset(data, T == T_value)
    d$n_surv <- d$y
    d$n_dead <- d$n - d$y
    fit <- tryCatch(
      stats::glm(cbind(n_surv, n_dead) ~ log10_t,
                 data = d,
                 family = stats::binomial("logit")),
      warning = function(w) w,
      error   = function(e) e
    )
    if (inherits(fit, "error") || inherits(fit, "warning") ||
        !inherits(fit, "glm") ||
        any(!is.finite(stats::coef(fit)))) {
      return(list(T = T_value, log10_LT50 = NA_real_,
                  se_log10_LT50 = NA_real_, success = FALSE))
    }
    co <- stats::coef(fit)
    if (co[["log10_t"]] >= 0 || !is.finite(co[["log10_t"]]) ||
        !is.finite(co[["(Intercept)"]])) {
      return(list(T = T_value, log10_LT50 = NA_real_,
                  se_log10_LT50 = NA_real_, success = FALSE))
    }
    # log10_LT50 = -b0 / b1   (where logit(p) = b0 + b1 * log10_t and p = 0.5)
    b0 <- co[["(Intercept)"]]; b1 <- co[["log10_t"]]
    log10_LT50 <- -b0 / b1
    V <- stats::vcov(fit)
    g <- c(-1 / b1, b0 / b1^2)                       # d(log10_LT50)/d(b0,b1)
    se_log10_LT50 <- sqrt(as.numeric(t(g) %*% V %*% g))
    list(T = T_value, log10_LT50 = log10_LT50,
         se_log10_LT50 = se_log10_LT50, success = TRUE)
  })

  s1 <- do.call(rbind.data.frame, per_temp)
  if (!all(s1$success) || sum(is.finite(s1$log10_LT50)) < 3L) {
    return(list(success = FALSE, stage1 = s1,
                z = list(point = NA_real_, lower = NA_real_, upper = NA_real_),
                CTmax_1hr = list(point = NA_real_, lower = NA_real_,
                                 upper = NA_real_)))
  }

  # Stage 2: OLS on (T, log10_LT50). Coefficients beta0_s2 (intercept on
  # natural T scale) and beta1_s2 (slope, °C^-1).
  s2 <- stats::lm(log10_LT50 ~ T, data = s1)
  co <- stats::coef(s2)
  V2 <- stats::vcov(s2)
  b0 <- co[["(Intercept)"]]; b1 <- co[["T"]]

  if (!is.finite(b1) || b1 >= 0) {
    return(list(success = FALSE, stage1 = s1,
                z = list(point = NA_real_, lower = NA_real_, upper = NA_real_),
                CTmax_1hr = list(point = NA_real_, lower = NA_real_,
                                 upper = NA_real_)))
  }

  # z = -1 / b1.  d(z)/d(b1) = 1 / b1^2.  SE: |d/dbeta| sqrt(V[b1, b1]).
  z_point <- -1 / b1
  z_se    <- abs(1 / b1^2) * sqrt(V2["T", "T"])

  # CTmax_1hr = (log10(t_ref_min) - b0) / b1
  log10_tref <- log10(t_ref_min)
  CTmax_pt   <- (log10_tref - b0) / b1
  g          <- c(-1 / b1, -(log10_tref - b0) / b1^2)  # d/d(b0, b1)
  CTmax_se   <- sqrt(as.numeric(t(g) %*% V2 %*% g))

  z_thr     <- stats::qnorm(0.975)
  list(
    success = TRUE,
    stage1  = s1,
    z = list(point = z_point,
             lower = z_point - z_thr * z_se,
             upper = z_point + z_thr * z_se,
             se    = z_se),
    CTmax_1hr = list(point = CTmax_pt,
                     lower = CTmax_pt - z_thr * CTmax_se,
                     upper = CTmax_pt + z_thr * CTmax_se,
                     se    = CTmax_se)
  )
}

# ----- joint-4PL wrapper for the simulation (no caching) ----------------------

#' Joint Bayesian 4PL fit + extract for one simulated dataset
#'
#' Thin wrapper around [standardize_data()] + [fit_4pl()] + [extract_tdt()]
#' configured for simulation speed: 2 chains, 2000 iterations (warmup 1000),
#' minimal console output, no on-disk caching. Returns the same three-piece
#' summary shape as [fit_two_stage_classical()] so the two estimators can be
#' compared row-for-row, **plus** the full per-draw posterior of $z$,
#' $CT_{max_{1hr}}$, and $T_{crit}$ (length `ndraws` each). The per-draw
#' vectors let downstream code pool posteriors across simulations or recompute
#' alternative summaries without re-fitting.
#'
#' @param data Output of [sim_twostage_dataset()].
#' @param ndraws Posterior draws to keep from `extract_tdt()`. Default 1000.
#' @param chains Number of MCMC chains. Default 2.
#' @param iter   Total iterations per chain (half are warmup). Default 2000.
#' @param seed   Sampler seed.
#' @return Named list with `$success`, `$z`, `$CTmax_1hr`, `$T_crit`,
#'         `$draws` (per-draw posterior of all three quantities, as a
#'         data frame with columns `.draw`, `z`, `CTmax_1hr`, `T_crit`),
#'         and `$diagnostics`.
#' @export
fit_joint_4pl_sim <- function(data,
                              ndraws = 1000,
                              chains = 2,
                              iter   = 2000,
                              seed   = 1L) {
  std <- standardize_data(data,
                          temp = "T", duration = "t",
                          n_total = "n", n_surv = "y",
                          duration_unit = "minutes")
  wf <- tryCatch(
    fit_4pl(std,
            chains     = chains,
            iter       = iter,
            cores      = chains,
            seed       = seed,
            refresh    = 0,
            silent     = 2,
            backend    = "cmdstanr",
            control    = list(adapt_delta = 0.9, max_treedepth = 12)),
    error = function(e) e
  )
  if (inherits(wf, "error")) {
    return(list(success = FALSE, error = conditionMessage(wf),
                z = list(point = NA_real_, lower = NA_real_, upper = NA_real_),
                CTmax_1hr = list(point = NA_real_, lower = NA_real_,
                                 upper = NA_real_),
                T_crit    = list(point = NA_real_, lower = NA_real_,
                                 upper = NA_real_),
                draws = NULL))
  }

  # The bias simulation generates beta-binomial *lethal* data, so opt into
  # the rate-multiplier T_crit for the joint-4PL summary.
  et <- tryCatch(
    suppressMessages(extract_tdt(wf, t_ref = 60, time_multiplier = 1,
                                  ndraws = ndraws, lethal = TRUE)),
    error = function(e) e
  )
  if (inherits(et, "error")) {
    return(list(success = FALSE, error = conditionMessage(et),
                z = list(point = NA_real_, lower = NA_real_, upper = NA_real_),
                CTmax_1hr = list(point = NA_real_, lower = NA_real_,
                                 upper = NA_real_),
                T_crit    = list(point = NA_real_, lower = NA_real_,
                                 upper = NA_real_),
                draws = NULL))
  }

  # Per-draw vectors paired by .draw index. The three extract_tdt slots all
  # carry .draw; inner_join keeps only draws that survived every primitive
  # (e.g. drops any draw whose LT50 OLS regression failed).
  z_d <- et$z$draws       |> dplyr::select(.draw, z)
  c_d <- et$CTmax$draws   |> dplyr::transmute(.draw, CTmax_1hr = temp)
  t_d <- et$T_crit$draws  |> dplyr::transmute(.draw, T_crit    = temp)
  draws_df <- z_d |>
    dplyr::inner_join(c_d, by = ".draw") |>
    dplyr::inner_join(t_d, by = ".draw")

  diag_tbl <- tryCatch(diagnose_tdt_fit(wf), error = function(e) NULL)
  rhat_max <- if (is.null(diag_tbl)) NA_real_ else diag_tbl$rhat_max
  divergences <- if (is.null(diag_tbl)) NA_integer_ else diag_tbl$divergences

  list(
    success = TRUE,
    z = list(point = et$z$summary$z_median,
             lower = et$z$summary$z_lower,
             upper = et$z$summary$z_upper),
    CTmax_1hr = list(point = et$CTmax$summary$temp_median,
                     lower = et$CTmax$summary$temp_lower,
                     upper = et$CTmax$summary$temp_upper),
    T_crit    = list(point = et$T_crit$summary$temp_median,
                     lower = et$T_crit$summary$temp_lower,
                     upper = et$T_crit$summary$temp_upper),
    draws = draws_df,
    diagnostics = list(rhat_max = rhat_max, divergences = divergences)
  )
}

# ----- per-simulation result extractor ----------------------------------------

#' Tidy one (joint 4PL, two-stage) pair against the truth
#'
#' Computes signed bias, interval coverage, and interval width for both
#' estimators and both quantities ($z$, $CT_{max_{1hr}}$). Returns a long-format
#' tibble — one row per (method, quantity).
#'
#' @param joint Output of [fit_joint_4pl_sim()].
#' @param ts    Output of [fit_two_stage_classical()].
#' @param truth Output of [sim_twostage_truth()].
#' @param sim_id Integer simulation index.
#' @param scenario Character scenario label (e.g. `"n3"`, `"n5"`).
#' @param runtime_sec Numeric — wall time spent on the joint 4PL fit.
#' @return Tibble with columns `sim_id`, `scenario`, `method`, `quantity`,
#'         `truth`, `estimate`, `bias`, `lower`, `upper`, `covered`,
#'         `width`, `success`.
#' @export
sim_twostage_result_row <- function(joint, ts, truth, sim_id, scenario,
                                    runtime_sec = NA_real_) {
  pack <- function(method, q_name, est, lo, hi, true_val, success) {
    tibble::tibble(
      sim_id   = sim_id,
      scenario = scenario,
      method   = method,
      quantity = q_name,
      truth    = true_val,
      estimate = est,
      bias     = est - true_val,
      lower    = lo,
      upper    = hi,
      covered  = is.finite(lo) & is.finite(hi) & lo <= true_val & true_val <= hi,
      width    = hi - lo,
      success  = success
    )
  }

  rows <- dplyr::bind_rows(
    pack("joint_4pl", "z",
         joint$z$point, joint$z$lower, joint$z$upper,
         truth$z_true, joint$success),
    pack("joint_4pl", "CTmax_1hr",
         joint$CTmax_1hr$point, joint$CTmax_1hr$lower, joint$CTmax_1hr$upper,
         truth$CTmax_1hr_true, joint$success),
    pack("two_stage", "z",
         ts$z$point, ts$z$lower, ts$z$upper,
         truth$z_true, ts$success),
    pack("two_stage", "CTmax_1hr",
         ts$CTmax_1hr$point, ts$CTmax_1hr$lower, ts$CTmax_1hr$upper,
         truth$CTmax_1hr_true, ts$success)
  )
  rows$runtime_sec <- runtime_sec
  rows
}
