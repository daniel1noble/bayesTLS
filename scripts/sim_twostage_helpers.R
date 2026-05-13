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
# Beta-binomial 4PL with temperature only on mid. Derived true z and CTmax_1hr
# are deterministic from these constants.

#' Default truth parameters for the two-stage bias simulation
#'
#' Beta-binomial 4PL with temperature varying only on the midpoint. Returns the
#' fixed truth used by [sim_twostage_dataset()] and by the analytical
#' z / CTmax_1hr targets the simulation judges bias against.
#'
#' @return Named list of truth parameters.
#' @export
sim_twostage_truth <- function() {
  out <- list(
    ell     = 0.05,
    u       = 0.92,
    k       = 8,
    m_beta0 = 1.5,
    m_beta1 = -0.15,
    T_bar   = 34,
    phi     = 5
  )
  out$z_true         <- -1 / out$m_beta1
  out$CTmax_1hr_true <- (
    log10(60) - out$m_beta0 -
      (1 / out$k) * log((out$u - 0.5) / (0.5 - out$ell))
  ) / out$m_beta1 + out$T_bar
  out
}

#' Default factorial design grid for the two-stage bias simulation
#'
#' Five assay temperatures and six durations matching the supplement's
#' simulation. Returns the (temp, duration) crossing only — replication is
#' added by [sim_twostage_dataset()] per scenario.
#'
#' @return Tibble with columns `T`, `t`, `log10_t`, `T_c`.
#' @export
sim_twostage_grid <- function() {
  truth <- sim_twostage_truth()
  tidyr::expand_grid(
    T = c(30, 32, 34, 36, 38),
    t = c(1, 5, 15, 45, 135, 405)
  ) |>
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
  grid <- sim_twostage_grid()
  design <- tidyr::expand_grid(grid, rep = seq_len(n_reps))

  p_true <- with(truth,
    ell + (u - ell) /
      (1 + exp(k * (design$log10_t - (m_beta0 + m_beta1 * design$T_c))))
  )
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

  et <- tryCatch(
    extract_tdt(wf, t_ref = 60, time_multiplier = 1, ndraws = ndraws),
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
