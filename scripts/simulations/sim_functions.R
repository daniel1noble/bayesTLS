# ============================================================================
# Building blocks for the two-stage-bias simulation.
#
# These are the reusable wrappers the runner (run_simulations.R) assembles.
# They are NOT part of the user-facing bayesTLS package — they are simulation
# internals private to this manuscript. They depend on bayesTLS's public API
# (standardize_data, fit_4pl, diagnose_tdt_fit; extract_joint() below replaces the
# retired extract_tdt()).
#
# The runner sources this file and then calls, per simulated dataset:
#   sim_dataset() -> fit_two_stage() x2 -> fit_joint_4pl() -> score_run()
# and aggregates the per-sim files with collect_raw() + the summary helpers.
#
# Nothing here parses CLI args or shells out — every step is plain R.
# ============================================================================

# ----------------------------------------------------------------------------
# 1. Study design: which (temperature x duration) cells a design measures
# ----------------------------------------------------------------------------

#' Temperature and duration grid for a named design
#'
#' A single lookup so the truth calculation and the data generator always agree
#' on the cells a design measures. All designs share the standard five assay
#' temperatures except `"sparse"`. The `tmax*` designs hold those temperatures
#' fixed and only shrink the *longest* exposure — they exist to ask how well a
#' shortened time window recovers z and CTmax_1hr (the design simulation).
#'
#' @param design One of `"full"`, `"sparse"`, `"tmax060"`, `"tmax120"`,
#'   `"tmax240"`, `"tmax405"`.
#' @return List with `temps` and `durations` (minutes).
design_spec <- function(design) {
  std_temps <- c(30, 32, 34, 36, 38)
  switch(design,
    full    = list(temps = std_temps,     durations = c(1, 5, 15, 45, 135, 405)),
    sparse  = list(temps = c(30, 34, 38), durations = c(1, 15, 135, 405)),
    # Max-exposure-time sweep: 6 log-spaced durations capped at 60/120/240/405
    # min. tmax405 reproduces the "full" grid (the good long-window design).
    tmax060 = list(temps = std_temps,     durations = c(1, 3, 8, 20, 40, 60)),
    tmax120 = list(temps = std_temps,     durations = c(1, 4, 12, 30, 70, 120)),
    tmax240 = list(temps = std_temps,     durations = c(1, 5, 18, 50, 130, 240)),
    tmax405 = list(temps = std_temps,     durations = c(1, 5, 15, 45, 135, 405)),
    # Extreme-temperature designs (2026-10-29 design question, see
    # notes/twostage_sim/): concentrate cells at the ends of the temperature
    # range instead of filling the middle. Durations span both ends' LT50
    # windows (t50 at 30 C ~ 126 min; at 38 C ~ 8 min under the baseline
    # truth); `extmid` keeps the centre temperature as a linearity check.
    ext2    = list(temps = c(30, 38),     durations = c(4, 10, 30, 60, 130, 300)),
    extmid  = list(temps = c(30, 34, 38), durations = c(4, 10, 30, 60, 130, 300)),
    stop("Unknown design: '", design, "'. Use full, sparse, tmax060/120/240/405, ",
         "ext2, or extmid.", call. = FALSE)
  )
}

#' Factorial (temperature x duration) grid for a design
#'
#' Returns the crossing only — replication is added by [sim_dataset()].
#'
#' @param design Design label (see [design_spec()]).
#' @return Tibble with columns `T`, `t`, `log10_t`.
sim_design_grid <- function(design = "full") {
  sp <- design_spec(design)
  tidyr::expand_grid(T = sp$temps, t = sp$durations) |>
    dplyr::mutate(log10_t = log10(t))
}

# ----------------------------------------------------------------------------
# 2. Truth: the data-generating parameters and the quantities to recover
# ----------------------------------------------------------------------------

#' Data-generating truth for the two-stage bias simulation
#'
#' Beta-binomial (or binomial) 4PL. Temperature always acts on the midpoint
#' (slope `m_beta1`); DGP variants add linear temperature slopes to the upper
#' asymptote, lower asymptote, and/or steepness:
#' \itemize{
#'   \item `"baseline"`  — only the midpoint varies with T.
#'   \item `"sym_ul"`    — u and ell shift symmetrically; (u+ell)/2 preserved.
#'   \item `"asym_u"`    — u shifts with T, ell constant.
#'   \item `"varying_k"` — k shifts with T; u, ell constant.
#' }
#'
#' The z and CTmax_1hr *targets* both estimators are scored against are the OLS
#' slope/intercept of log10(LT50_{p=0.5})(T) at the design temperatures, taken
#' from the (possibly T-varying) truth surface — see [compute_ols_truth()].
#'
#' @param dgp One of `"baseline"`, `"sym_ul"`, `"asym_u"`, `"varying_k"`.
#' @param u_0,ell_0 Optional overrides for the asymptote intercepts (value at
#'   T = T_bar). `NULL` keeps the preset.
#' @param u_beta1,ell_beta1,k_beta1 Optional overrides for the temperature
#'   slopes of u, ell, k. `NULL` keeps the preset.
#' @param batch_sd SD of a random shift of each day's midpoints (log10 min), drawn anew for
#'   every dataset. 0 (the default) gives the paper's main scenarios.
#' @param schedule Which cups a day holds (see `fig_batch_sequence.R`): `"A"` one temperature
#'   per day, `"B"` every temperature every day, with replicate cup r on day r, `"C"` each
#'   temperature on two days (one temperature a day; its first ceiling(n_reps / 2) replicate cups
#'   on one day and the rest on another).
#' @param design Design label; controls the grid the OLS truth is evaluated on.
#' @param family `"beta_binomial"` (default, overdispersed) or `"binomial"`.
#'   Travels with the truth so a forgotten override can't silently switch
#'   likelihoods between cells.
#' @return Named list of truth parameters plus `z_true`, `CTmax_1hr_true`.
sim_truth <- function(dgp       = "baseline",
                      u_0       = NULL,
                      ell_0     = NULL,
                      u_beta1   = NULL,
                      ell_beta1 = NULL,
                      k_beta1   = NULL,
                      batch_sd  = 0,
                      schedule  = c("A", "B", "C"),
                      design    = "full",
                      family    = c("beta_binomial", "binomial")) {
  family <- match.arg(family)
  base <- list(
    ell = 0.05, u = 0.92, k = 8,
    m_beta0 = 1.5, m_beta1 = -0.15, T_bar = 34, phi = 5,
    u_beta1 = 0, ell_beta1 = 0, k_beta1 = 0
  )
  out <- switch(dgp,
    baseline  = base,
    sym_ul    = utils::modifyList(base, list(u_beta1 = -0.01, ell_beta1 = 0.01)),
    asym_u    = utils::modifyList(base, list(u_beta1 = -0.01)),
    varying_k = utils::modifyList(base, list(k_beta1 = 0.25)),
    stop("Unknown dgp: '", dgp, "'. ",
         "Use baseline, sym_ul, asym_u, or varying_k.", call. = FALSE)
  )
  if (!is.null(u_0))       out$u         <- u_0
  if (!is.null(ell_0))     out$ell       <- ell_0
  if (!is.null(u_beta1))   out$u_beta1   <- u_beta1
  if (!is.null(ell_beta1)) out$ell_beta1 <- ell_beta1
  if (!is.null(k_beta1))   out$k_beta1   <- k_beta1
  # Day (batch) effects for the non-independence simulation (2026-11-28; off in the paper's
  # main scenarios, whose truth lists stay unchanged).
  schedule <- match.arg(schedule)
  if (batch_sd > 0)    out$batch_sd <- batch_sd
  if (schedule != "A") out$schedule <- schedule

  out$dgp    <- dgp
  out$design <- design
  out$family <- family

  tt <- compute_ols_truth(out, design = design)
  out$z_true         <- tt$z_true
  out$CTmax_1hr_true <- tt$CTmax_1hr_true
  out
}

#' OLS z / CTmax_1hr targets at the design temperatures
#'
#' For each design temperature, evaluate
#' `log10(LT50_{p=0.5})(T) = mid(T) + (1/k(T)) * log((u(T)-0.5)/(0.5-ell(T)))`
#' on the (possibly T-varying) truth surface, then OLS-fit a line and read off
#' `z_true = -1/slope` and `CTmax_1hr_true = (log10(60) - intercept)/slope`.
#' The duration window does not enter here, so the truth is identical across the
#' tmax sweep — only a design's ability to recover it changes.
#'
#' @param p Truth list (from [sim_truth()] before z/CTmax are written).
#' @param design Design label (supplies the temperatures).
#' @return List with `z_true` and `CTmax_1hr_true`.
compute_ols_truth <- function(p, design = "full") {
  Ts    <- design_spec(design)$temps
  T_c   <- Ts - p$T_bar
  u_T   <- p$u   + p$u_beta1   * T_c
  ell_T <- p$ell + p$ell_beta1 * T_c
  k_T   <- p$k   + p$k_beta1   * T_c
  mid_T <- p$m_beta0 + p$m_beta1 * T_c

  # Feasibility: need 0 < ell < 0.5 < u < 1 and k > 0 at every design temp so
  # log10(LT50_{p=0.5}) is defined and the beta-binomial shapes are positive.
  if (any(u_T <= 0.5) || any(ell_T >= 0.5) || any(k_T <= 0) ||
      any(u_T >= 1)   || any(ell_T <= 0))
    stop("DGP gives infeasible u/ell/k at design temperatures: ",
         "need 0 < ell < 0.5 < u < 1 and k > 0. u_T = ",
         paste(sprintf("%.3f", u_T), collapse = ", "), "; ell_T = ",
         paste(sprintf("%.3f", ell_T), collapse = ", "), "; k_T = ",
         paste(sprintf("%.3f", k_T), collapse = ", "), ".", call. = FALSE)

  log10_lt50 <- mid_T + (1 / k_T) * log((u_T - 0.5) / (0.5 - ell_T))
  co    <- stats::coef(stats::lm(log10_lt50 ~ Ts))
  slope <- co[["Ts"]]; intercept <- co[["(Intercept)"]]
  list(z_true         = -1 / slope,
       CTmax_1hr_true = (log10(60) - intercept) / slope)
}

# ----------------------------------------------------------------------------
# 3. Data generator
# ----------------------------------------------------------------------------

#' Simulate one 4PL dataset (binomial or beta-binomial)
#'
#' `n_reps` replicate cups per (temperature x duration) cell, each with N
#' individuals drawn uniformly from `n_ind_range`. Under `"beta_binomial"` the
#' per-cup probability is drawn from Beta(p*phi, (1-p)*phi) then the count from
#' Binomial(N, p_cup) — capturing cup-to-cup overdispersion. Under `"binomial"`
#' the count is Binomial(N, p) directly.
#'
#' `truth` is REQUIRED (no default): a forgotten override must not silently fall
#' back to the baseline DGP (the 2026-05-15 incident in feedback_sim_preflight).
#'
#' @param n_reps Replicate cups per cell.
#' @param seed Integer seed.
#' @param truth Truth list from [sim_truth()].
#' @param n_ind_range Length-2 range for N per cup. Default `c(10, 20)`.
#' @return Tibble: `T`, `t`, `log10_t`, `T_c`, `rep`, `day`, `shift` (the cup's day shift, log10 min),
#'   `n`, `y`, `p_true`.
sim_dataset <- function(n_reps, seed, truth, n_ind_range = c(10, 20)) {
  if (missing(truth))
    stop("sim_dataset(): `truth` is required; pass sim_truth(...).", call. = FALSE)
  set.seed(seed)
  grid   <- sim_design_grid(design = truth$design %||% "full")
  design <- tidyr::expand_grid(grid, rep = seq_len(n_reps))

  design$T_c <- design$T - truth$T_bar
  u_T   <- truth$u   + truth$u_beta1   * design$T_c
  ell_T <- truth$ell + truth$ell_beta1 * design$T_c
  k_T   <- truth$k   + truth$k_beta1   * design$T_c
  mid_T <- truth$m_beta0 + truth$m_beta1 * design$T_c
  # Day of each cup. Each day shifts the midpoint of every cup it holds.
  Ts <- sort(unique(design$T))
  design$day <- switch(truth$schedule %||% "A",
    A = match(design$T, Ts),                                     # one temperature per day
    B = design$rep,                                              # every temperature every day
    C = match(design$T, Ts) + length(Ts) * (design$rep > ceiling(n_reps / 2)))  # each temperature on two days
  design$shift <- if ((truth$batch_sd %||% 0) > 0) stats::rnorm(max(design$day), 0, truth$batch_sd)[design$day] else 0
  mid_T <- mid_T + design$shift
  p_true <- ell_T + (u_T - ell_T) / (1 + exp(k_T * (design$log10_t - mid_T)))

  n   <- sample(seq(n_ind_range[1], n_ind_range[2]), nrow(design), replace = TRUE)
  fam <- truth$family %||% "beta_binomial"
  y <- if (fam == "beta_binomial") {
    p_cup <- stats::rbeta(nrow(design), p_true * truth$phi, (1 - p_true) * truth$phi)
    stats::rbinom(nrow(design), size = n, prob = p_cup)
  } else if (fam == "binomial") {
    stats::rbinom(nrow(design), size = n, prob = p_true)
  } else {
    stop("Unknown family '", fam, "'.", call. = FALSE)
  }

  design$n <- n; design$y <- y; design$p_true <- p_true
  tibble::as_tibble(design)
}

# ----------------------------------------------------------------------------
# 4. Estimators
# ----------------------------------------------------------------------------

#' Classical two-stage TDT estimator (per-temperature GLM + OLS)
#'
#' Stage 1: at each temperature fit a dose-response GLM of survival on
#' `log10(duration)` and read off `log10(LT50)` at p = 0.5. Stage 2: OLS of
#' `log10(LT50)` on temperature, giving `z = -1/beta1` and
#' `CTmax_1hr = (log10(60) - beta0)/beta1`. 95% CIs come from the delta method
#' on the Stage-2 coefficients, returned with both Normal- and t-quantiles
#' (the t-quantile is the small-sample fix for the few Stage-2 residual df).
#'
#' The `stage1` argument selects the Stage-1 likelihood — this is the single
#' knob that used to be two near-identical functions:
#' \itemize{
#'   \item `"binomial"`     — `stats::glm` logit binomial (the field default).
#'   \item `"betabinomial"` — `glmmTMB::glmmTMB` beta-binomial (overdispersion
#'         handled at Stage 1; isolates the likelihood cost from the
#'         architecture cost).
#' }
#'
#' @param data Output of [sim_dataset()].
#' @param stage1 `"binomial"` or `"betabinomial"`.
#' @param t_ref_min Reference exposure time (minutes). Default 60.
#' @return List: `success`, `stage1` table, `z` and `CTmax_1hr` each as
#'   `list(point, lower, upper, lower_t, upper_t, se)`.
fit_two_stage <- function(data, stage1 = c("binomial", "betabinomial"),
                          t_ref_min = 60) {
  stage1 <- match.arg(stage1)
  if (stage1 == "betabinomial" && !requireNamespace("glmmTMB", quietly = TRUE))
    stop("fit_two_stage(stage1 = 'betabinomial') needs the glmmTMB package.",
         call. = FALSE)

  fail <- list(success = FALSE,
               z         = list(point = NA_real_, lower = NA_real_, upper = NA_real_),
               CTmax_1hr = list(point = NA_real_, lower = NA_real_, upper = NA_real_))

  # ---- Stage 1: one GLM per temperature; read off log10(LT50) + its SE ----
  per_temp <- lapply(sort(unique(data$T)), function(T_value) {
    d <- subset(data, T == T_value)
    d$n_surv <- d$y; d$n_dead <- d$n - d$y
    bad <- list(T = T_value, log10_LT50 = NA_real_, se_log10_LT50 = NA_real_,
                success = FALSE)

    if (stage1 == "binomial") {
      # Benign "fitted probabilities 0/1" warnings near asymptotes are expected;
      # genuine failures are caught by the finiteness / negative-slope checks.
      fit <- tryCatch(suppressWarnings(
        stats::glm(cbind(n_surv, n_dead) ~ log10_t, data = d,
                   family = stats::binomial("logit"))), error = function(e) e)
      if (inherits(fit, "error") || !inherits(fit, "glm") ||
          any(!is.finite(stats::coef(fit)))) return(bad)
      co <- stats::coef(fit); V <- stats::vcov(fit)
    } else {
      fit <- tryCatch(suppressWarnings(suppressMessages(
        glmmTMB::glmmTMB(cbind(n_surv, n_dead) ~ log10_t, data = d,
                         family = glmmTMB::betabinomial(link = "logit")))),
        error = function(e) e)
      if (inherits(fit, "error") || !inherits(fit, "glmmTMB")) return(bad)
      co <- tryCatch(glmmTMB::fixef(fit)$cond, error = function(e) NULL)
      V  <- tryCatch(stats::vcov(fit)$cond,    error = function(e) NULL)
      if (is.null(co) || is.null(V) || any(!is.finite(co))) return(bad)
    }

    if (co[["log10_t"]] >= 0 || !is.finite(co[["log10_t"]]) ||
        !is.finite(co[["(Intercept)"]])) return(bad)
    b0 <- co[["(Intercept)"]]; b1 <- co[["log10_t"]]
    log10_LT50 <- -b0 / b1
    g  <- c(-1 / b1, b0 / b1^2)                       # d(log10_LT50)/d(b0, b1)
    list(T = T_value, log10_LT50 = log10_LT50,
         se_log10_LT50 = sqrt(as.numeric(t(g) %*% V %*% g)), success = TRUE)
  })

  s1 <- do.call(rbind.data.frame, per_temp)
  if (!all(s1$success) || sum(is.finite(s1$log10_LT50)) < 3L)
    return(utils::modifyList(fail, list(stage1 = s1)))

  # ---- Stage 2: unweighted OLS of log10(LT50) on temperature ----
  s2 <- stats::lm(log10_LT50 ~ T, data = s1)
  co <- stats::coef(s2); V2 <- stats::vcov(s2)
  b0 <- co[["(Intercept)"]]; b1 <- co[["T"]]
  if (!is.finite(b1) || b1 >= 0)
    return(utils::modifyList(fail, list(stage1 = s1)))

  z_point  <- -1 / b1
  z_se     <- abs(1 / b1^2) * sqrt(V2["T", "T"])
  log10_tr <- log10(t_ref_min)
  CTmax_pt <- (log10_tr - b0) / b1
  gC       <- c(-1 / b1, -(log10_tr - b0) / b1^2)     # d(CTmax)/d(b0, b1)
  CTmax_se <- sqrt(as.numeric(t(gC) %*% V2 %*% gC))

  # Stage-2 OLS has n_T - 2 residual df (3 for full, 1 for sparse). Normal
  # quantiles under-cover at small df; t-quantiles are the fix. Return both.
  df_resid <- stats::df.residual(s2)
  z_norm   <- stats::qnorm(0.975)
  z_t      <- if (is.finite(df_resid) && df_resid > 0)
                stats::qt(0.975, df = df_resid) else NA_real_

  ci <- function(pt, se) list(point = pt,
    lower = pt - z_norm * se, upper = pt + z_norm * se,
    lower_t = pt - z_t * se,  upper_t = pt + z_t * se, se = se)

  list(success = TRUE, stage1 = s1, df_resid = df_resid,
       z = ci(z_point, z_se), CTmax_1hr = ci(CTmax_pt, CTmax_se))
}

#' Joint 4PL extraction: z, CTmax_1hr and T_crit at one threshold
#'
#' Reproduces, draw for draw, the retired `bayesTLS::extract_tdt()` that
#' produced the raw results (the mid-June 2026 code, e.g. git
#' `977498a:R/extract_tdt.R`; its last version, `19c123f^:R/extract_tdt.R`,
#' gives the same results), including its use of the random-number stream: the
#' draw subsample, the subsample its LT-curve step drew next, and the T_crit
#' rate draws. Two changes (2026-11-26), both at the ABSOLUTE threshold only:
#' (1) extract_tdt() inverted each draw's LT curve on the fixed grid
#' `min(T) - 2` to `max(T) + 2` C and silently dropped every draw whose
#' CTmax_1hr crossing lay outside it, so the summary came from a truncated
#' posterior; here such a draw is searched for beyond the grid as bayesTLS's
#' tls() does, the first crossing outward on the same 0.05 C spacing (`reach`
#' C each side). (2) A draw whose LT curve is not monotone over the grid now
#' takes its single falling crossing, or none when it has none or several
#' (extract_tdt() returned a y-sorted `approx()` value, which need not be a
#' crossing at all), and a target outside a draw's fitted asymptotes is
#' treated as an LT of zero or infinity, with the crossing placed on the
#' survival curve. Draws still without a single falling crossing within
#' `reach` (for example, an LT that lengthens with temperature) are left out
#' of the summary and counted in `n_uncrossed`. Monotone draws crossing on the
#' grid, the absolute z and the relative threshold are unchanged.
#'
#' @param wf Fitted `bayes_tls` workflow (single condition, minutes, t_ref 60).
#' @param target_surv `"relative"` or `"absolute"` (literal 50% survival).
#' @param ndraws Posterior draws to subsample (the random stream is used
#'   exactly as extract_tdt() used it).
#' @param TC_rate_range HI-rate range for the T_crit rate draws.
#' @param reach Degrees C searched beyond each end of the grid (absolute only).
#' @return List with `z`, `CTmax`, `T_crit` (each `draws` + `summary`, the
#'   extract_tdt() shapes), `n_uncrossed` (absolute draws still without a
#'   crossing; 0 for relative) and `draw_ids` (the posterior draws used; NULL =
#'   all).
extract_joint <- function(wf, target_surv = c("relative", "absolute"),
                          ndraws = 1000, TC_rate_range = c(0.1, 1), reach = 30) {
  target_surv <- match.arg(target_surv)
  ns     <- asNamespace("bayesTLS")
  fit    <- get_brmsfit(wf)
  data   <- wf$data
  Tbar   <- wf$meta$temp_mean
  bnd    <- wf$meta$bounds
  abs_p  <- identical(target_surv, "absolute")
  p      <- if (abs_p) 0.5 else NA_real_
  target <- log10(60)                       # t_ref = 60 min; the model is in minutes
  trange <- range(data$temp, na.rm = TRUE)
  temp_grid <- seq(trange[1] - 2, trange[2] + 2, by = 0.05)    # extract_tdt()'s grid
  z_grid <- sort(unique(data$temp)); z_grid <- z_grid[is.finite(z_grid)]
  Lz <- length(z_grid); Lc <- length(temp_grid); h <- 1e-3

  # Random stream, in extract_tdt()'s order: (1) the draw subsample ...
  did <- ns$tls_draw_ids(fit, ndraws)
  combo <- c(0, z_grid - Tbar - h, z_grid - Tbar + h, temp_grid - Tbar)
  nd <- ns$tls_build_grid(fit$data, temp = "temp_c", temp_grid = combo)
  sp <- ns$tls_eval_subpars(fit, nd, bnd, draw_ids = did, mode = target_surv, p = p)
  # ... (2) the subsample its LT-curve step (derive_tdt_curve(ndraws = ndraws))
  # drew through brms::posterior_linpred(), i.e. sample(seq_len(n), ndraws) ...
  n_all <- brms::ndraws(fit)
  invisible(sample.int(n_all, min(ndraws, n_all)))

  c_minus <- 1L + seq_len(Lz); c_plus <- 1L + Lz + seq_len(Lz)
  c_ct    <- 1L + 2L * Lz + seq_len(Lc)
  zc <- ns$tls_zct_draws(sp$logLT[, c_plus, drop = FALSE], sp$logLT[, c_minus, drop = FALSE],
                         h, z_grid, sp$logLT[, c_ct, drop = FALSE], temp_grid,
                         sp$logLT[, 1], target = target, target_1hr = target, Tbar = Tbar,
                         use_closed = !abs_p, want_ctmax = !abs_p, want_ct1 = !abs_p)
  n_uncrossed <- 0L
  if (abs_p) {
    model <- ns$tdt_model(fit)
    spc <- lapply(sp[c("low", "up", "k", "mid", "logLT")], function(m) m[, c_ct, drop = FALSE])
    L  <- ns$tls_clamp_logLT(spc, p)
    Tc <- ns$tls_invert_logLT(L, target, temp_grid, R = ns$tls_surv_resid(spc, target, p, model))
    need <- which(!is.finite(Tc) &
                  vapply(seq_len(nrow(L)), function(i) ns$tls_outward_side(L[i, ] - target) != 0L,
                         logical(1)))
    if (length(need)) {
      n_out <- round(reach / 0.05)
      lower <- temp_grid[1]  - 0.05 * rev(seq_len(n_out))
      upper <- temp_grid[Lc] + 0.05 * seq_len(n_out)
      ndx <- nd[rep(1L, 2L * n_out), , drop = FALSE]
      ndx$temp_c <- c(lower, upper) - Tbar
      rownames(ndx) <- NULL
      spe <- ns$tls_eval_subpars(fit, ndx, bnd, draw_ids = if (is.null(did)) need else did[need],
                                 mode = "absolute", p = p, warn = FALSE)
      spf <- lapply(c(low = "low", up = "up", k = "k", mid = "mid", logLT = "logLT"), function(nm)
        cbind(spe[[nm]][, seq_len(n_out), drop = FALSE], spc[[nm]][need, , drop = FALSE],
              spe[[nm]][, n_out + seq_len(n_out), drop = FALSE]))
      Tc[need] <- ns$tls_invert_outward(ns$tls_clamp_logLT(spf, p), target,
                                        c(lower, temp_grid, upper), n_out + seq_len(Lc),
                                        R = ns$tls_surv_resid(spf, target, p, model))
    }
    n_uncrossed <- sum(!is.finite(Tc))
    zc$ctmax <- Tc; zc$ct1 <- Tc            # t_ref = 60: the 1-hour anchor is CTmax_1hr
  }

  qs <- function(v) stats::quantile(v, c(0.025, 0.5, 0.975), na.rm = TRUE, names = FALSE)
  cd <- tibble::tibble(.draw = seq_along(zc$ctmax), temp = zc$ctmax) |>
    dplyr::filter(is.finite(temp))
  cq <- qs(cd$temp)
  # ... and (3) the T_crit rate draws, one per draw with a finite z and anchor.
  paired <- dplyr::inner_join(dplyr::select(zc$z$draws, .draw, z),
                              tibble::tibble(.draw = seq_along(zc$ct1), CTmax_1hr = zc$ct1) |>
                                dplyr::filter(is.finite(CTmax_1hr)), by = ".draw")
  paired$log10_rate <- stats::runif(nrow(paired), log10(TC_rate_range[1] / 100),
                                    log10(TC_rate_range[2] / 100))
  tcd <- tibble::tibble(.draw = paired$.draw, temp = paired$CTmax_1hr + paired$z * paired$log10_rate,
                        log10_rate = paired$log10_rate)
  tq <- qs(tcd$temp)
  list(z = zc$z,
       CTmax = list(draws = cd, summary = tibble::tibble(temp_lower = cq[1], temp_median = cq[2],
                                                         temp_upper = cq[3])),
       T_crit = list(draws = tcd, summary = tibble::tibble(temp_lower = tq[1], temp_median = tq[2],
                                                           temp_upper = tq[3])),
       n_uncrossed = n_uncrossed, draw_ids = did)
}

#' Joint Bayesian 4PL fit + extract for one simulated dataset
#'
#' Wraps [standardize_data()] + [fit_4pl()] + [extract_joint()] (both the relative
#' and absolute survival thresholds on the same fit). Returns the same
#' point/CI shape as [fit_two_stage()] for z and CTmax_1hr, plus the per-draw
#' posterior (z, CTmax_1hr, T_crit), sampler diagnostics, and the fit runtime.
#'
#' @param data Output of [sim_dataset()].
#' @param ndraws Posterior draws to keep. Default 1000.
#' @param chains,iter,warmup,max_treedepth,adapt_delta Sampler settings.
#' @param seed Sampler seed.
#' @return Named list (see body); `success = FALSE` carries the error message.
fit_joint_4pl <- function(data, ndraws = 1000,
                          chains = 3, iter = 3000, warmup = 1500,
                          max_treedepth = 16, adapt_delta = 0.95, seed = 1L) {
  empty <- list(point = NA_real_, lower = NA_real_, upper = NA_real_)
  fail  <- function(msg) list(
    success = FALSE, error = msg,
    z = empty, z_abs = empty, CTmax_1hr = empty, CTmax_1hr_abs = empty,
    T_crit = empty, draws = NULL, draws_abs = NULL,
    diagnostics = NULL, runtime_sec = NA_real_)

  std <- standardize_data(data, temp = "T", duration = "t",
                          n_total = "n", n_surv = "y", duration_unit = "minutes")

  # max_treedepth = 16 (above the fit_4pl default 12): boundary-asymptote truths
  # (e.g. u = 0.999) sit at the inv_logit reparameterised bound and need long
  # NUTS paths. 3 chains x (3000 - 1500) = 4500 post-warmup draws.
  t0 <- Sys.time()
  wf <- tryCatch(
    fit_4pl(std, chains = chains, iter = iter, warmup = warmup, cores = chains,
            seed = seed, refresh = 0, silent = 2, backend = "cmdstanr",
            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)),
    error = function(e) e)
  runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(wf, "error")) return(utils::modifyList(fail(conditionMessage(wf)),
                                                      list(runtime_sec = runtime_sec)))

  # Both threshold definitions on the same fit. "relative" = midpoint between
  # the estimated asymptotes (extract_tdt default); "absolute" = literal 50%
  # survival (matches the analytical OLS truth). Second extract is sub-second.
  et_rel <- tryCatch(suppressMessages(extract_joint(wf, "relative", ndraws = ndraws)),
                     error = function(e) e)
  et_abs <- tryCatch(suppressMessages(extract_joint(wf, "absolute", ndraws = ndraws)),
                     error = function(e) e)
  if (inherits(et_rel, "error") || inherits(et_abs, "error")) {
    msg <- if (inherits(et_rel, "error")) conditionMessage(et_rel)
           else conditionMessage(et_abs)
    return(utils::modifyList(fail(msg), list(runtime_sec = runtime_sec)))
  }

  pack_draws <- function(et)
    et$z$draws |> dplyr::select(.draw, z) |>
      dplyr::inner_join(et$CTmax$draws  |> dplyr::transmute(.draw, CTmax_1hr = temp), ".draw") |>
      dplyr::inner_join(et$T_crit$draws |> dplyr::transmute(.draw, T_crit = temp),    ".draw")

  diag_tbl <- tryCatch(diagnose_tdt_fit(wf), error = function(e) NULL)

  list(
    success = TRUE, runtime_sec = runtime_sec,
    z = list(point = et_rel$z$summary$z_median,
             lower = et_rel$z$summary$z_lower, upper = et_rel$z$summary$z_upper),
    z_abs = list(point = et_abs$z$summary$z_median,
                 lower = et_abs$z$summary$z_lower, upper = et_abs$z$summary$z_upper),
    CTmax_1hr = list(point = et_rel$CTmax$summary$temp_median,
                     lower = et_rel$CTmax$summary$temp_lower,
                     upper = et_rel$CTmax$summary$temp_upper),
    CTmax_1hr_abs = list(point = et_abs$CTmax$summary$temp_median,
                         lower = et_abs$CTmax$summary$temp_lower,
                         upper = et_abs$CTmax$summary$temp_upper),
    T_crit = list(point = et_rel$T_crit$summary$temp_median,
                  lower = et_rel$T_crit$summary$temp_lower,
                  upper = et_rel$T_crit$summary$temp_upper),
    draws = pack_draws(et_rel), draws_abs = pack_draws(et_abs),
    abs_uncrossed = et_abs$n_uncrossed,
    diagnostics = diag_tbl)
}

# ----------------------------------------------------------------------------
# 5. Scoring: bundle one simulation's per-sim outputs
# ----------------------------------------------------------------------------

# Default diagnostics row when a joint fit failed or produced none.
.diag_default <- function()
  tibble::tibble(rhat_max = NA_real_, ess_bulk_min = NA_real_,
                 ess_tail_min = NA_real_, divergences = NA_integer_,
                 treedepth_hits = NA_integer_, bfmi_min = NA_real_,
                 rhat_pass = NA, ess_pass = NA, divergence_pass = NA,
                 treedepth_pass = NA, bfmi_pass = NA, all_pass = NA)

#' Score one simulation against the truth
#'
#' Produces the four objects that make up a per-sim result file:
#' \itemize{
#'   \item `row` — long bias/coverage rows, one per (method, quantity), for
#'     z and CTmax_1hr across all six method labels.
#'   \item `meta` — sim id, seed, runtime, success flags, sampler diagnostics.
#'   \item `draws`, `draws_abs` — the per-draw joint-4PL posteriors, tagged.
#' }
#' Method labels (unchanged from the previous harness, so the manuscript figures
#' keep working): `joint_4pl`, `joint_4pl_abs`, `two_stage_bin`,
#' `two_stage_bin_t`, `two_stage_bb`, `two_stage_bb_t`.
#'
#' @param joint Output of [fit_joint_4pl()].
#' @param ts_bin,ts_bb Output of [fit_two_stage()] (binomial / betabinomial).
#' @param truth Output of [sim_truth()].
#' @param sim_id,scenario,seed Identifiers for this run.
#' @return List with `row`, `meta`, `draws`, `draws_abs`.
score_run <- function(joint, ts_bin, ts_bb, truth, sim_id, scenario, seed = NA_integer_) {
  pack <- function(method, q, est, lo, hi, true_val, success)
    tibble::tibble(sim_id = sim_id, scenario = scenario, method = method,
                   quantity = q, truth = true_val, estimate = est, bias = est - true_val,
                   lower = lo, upper = hi,
                   covered = is.finite(lo) & is.finite(hi) & lo <= true_val & true_val <= hi,
                   width = hi - lo, success = success)

  row <- dplyr::bind_rows(
    pack("joint_4pl", "z", joint$z$point, joint$z$lower, joint$z$upper,
         truth$z_true, joint$success),
    pack("joint_4pl", "CTmax_1hr", joint$CTmax_1hr$point, joint$CTmax_1hr$lower,
         joint$CTmax_1hr$upper, truth$CTmax_1hr_true, joint$success),
    pack("joint_4pl_abs", "z", joint$z_abs$point, joint$z_abs$lower,
         joint$z_abs$upper, truth$z_true, joint$success),
    pack("joint_4pl_abs", "CTmax_1hr", joint$CTmax_1hr_abs$point,
         joint$CTmax_1hr_abs$lower, joint$CTmax_1hr_abs$upper,
         truth$CTmax_1hr_true, joint$success),
    pack("two_stage_bin", "z", ts_bin$z$point, ts_bin$z$lower, ts_bin$z$upper,
         truth$z_true, ts_bin$success),
    pack("two_stage_bin", "CTmax_1hr", ts_bin$CTmax_1hr$point,
         ts_bin$CTmax_1hr$lower, ts_bin$CTmax_1hr$upper,
         truth$CTmax_1hr_true, ts_bin$success),
    pack("two_stage_bin_t", "z", ts_bin$z$point, ts_bin$z$lower_t,
         ts_bin$z$upper_t, truth$z_true, ts_bin$success),
    pack("two_stage_bin_t", "CTmax_1hr", ts_bin$CTmax_1hr$point,
         ts_bin$CTmax_1hr$lower_t, ts_bin$CTmax_1hr$upper_t,
         truth$CTmax_1hr_true, ts_bin$success),
    pack("two_stage_bb", "z", ts_bb$z$point, ts_bb$z$lower, ts_bb$z$upper,
         truth$z_true, ts_bb$success),
    pack("two_stage_bb", "CTmax_1hr", ts_bb$CTmax_1hr$point,
         ts_bb$CTmax_1hr$lower, ts_bb$CTmax_1hr$upper,
         truth$CTmax_1hr_true, ts_bb$success),
    pack("two_stage_bb_t", "z", ts_bb$z$point, ts_bb$z$lower_t,
         ts_bb$z$upper_t, truth$z_true, ts_bb$success),
    pack("two_stage_bb_t", "CTmax_1hr", ts_bb$CTmax_1hr$point,
         ts_bb$CTmax_1hr$lower_t, ts_bb$CTmax_1hr$upper_t,
         truth$CTmax_1hr_true, ts_bb$success)
  )
  row$runtime_sec <- joint$runtime_sec %||% NA_real_

  diag_row <- if (isTRUE(joint$success) && !is.null(joint$diagnostics))
                joint$diagnostics else .diag_default()
  meta <- tibble::tibble(sim_id = sim_id, scenario = scenario, seed = seed,
                         joint_sec = joint$runtime_sec %||% NA_real_,
                         joint_ok = isTRUE(joint$success),
                         ts_bin_ok = isTRUE(ts_bin$success),
                         ts_bb_ok = isTRUE(ts_bb$success),
                         abs_uncrossed = joint$abs_uncrossed %||% NA_integer_,
                         extraction = "extract_joint 2026-11-26") |>
    dplyr::bind_cols(diag_row)

  tag <- function(d) if (is.null(d)) NULL else
    dplyr::mutate(d, sim_id = sim_id, scenario = scenario) |>
      dplyr::select(sim_id, scenario, .draw, z, CTmax_1hr, T_crit)

  list(row = row, meta = meta,
       draws = if (joint$success) tag(joint$draws) else NULL,
       draws_abs = if (joint$success) tag(joint$draws_abs) else NULL)
}

# ----------------------------------------------------------------------------
# 6. Aggregation: per-sim files -> the 8 scenario output objects
# ----------------------------------------------------------------------------

#' One simulation, start to finish
#'
#' Simulate one dataset -> fit two-stage (both Stage-1 likelihoods) -> fit the
#' joint 4PL -> score all three against the truth -> save. Resume-safe: a sim
#' whose raw file exists is skipped unless `force`.
#'
#' @param sim_id Simulation id within the scenario.
#' @param sc One scenario row (label, index, dgp, family, n_reps, design, ...).
#' @param truth Output of [sim_truth()] for `sc`.
#' @param raw_dir Directory for the per-sim raw files.
#' @param sampler Sampler settings for [fit_joint_4pl()].
#' @param master_seed Seed base: seed = master_seed + 1000 * sc$index + sim_id.
#' @param force Refit even if the raw file exists.
run_one_sim <- function(sim_id, sc, truth, raw_dir, sampler, master_seed, force = FALSE) {
  f <- file.path(raw_dir, sprintf("sim_%05d.rds", sim_id))
  if (!force && file.exists(f)) return(invisible(NULL))
  seed   <- master_seed + 1000L * sc$index + sim_id
  data   <- sim_dataset(n_reps = sc$n_reps, seed = seed, truth = truth)
  ts_bin <- fit_two_stage(data, stage1 = "binomial")
  ts_bb  <- fit_two_stage(data, stage1 = "betabinomial")
  joint  <- do.call(fit_joint_4pl, c(list(data, seed = seed), sampler))
  save_retry(score_run(joint, ts_bin, ts_bb, truth, sim_id, sc$label, seed), f)
  invisible(NULL)
}

#' Keep the June two-stage results in a re-run simulation (2026-11-26)
#'
#' The re-run reproduces the joint 4PL exactly (same software) but the glmmTMB
#' installed now differs from June's on a few datasets, and the two-stage arm
#' has no bug. So the re-run result keeps the original file's two-stage rows
#' and success flags; its joint rows are new. Also returns a one-row check of
#' what reproduced (largest absolute difference in estimate or interval).
#'
#' @param new,old The re-run and the original per-sim results ([score_run()]).
#' @return List with `result` (merged per-sim result) and `check` (tibble).
merge_rerun <- function(new, old) {
  key  <- paste(new$row$method, new$row$quantity)
  m    <- match(key, paste(old$row$method, old$row$quantity))
  d    <- pmax(abs(new$row$estimate - old$row$estimate[m]),
               abs(new$row$lower - old$row$lower[m]),
               abs(new$row$upper - old$row$upper[m]))
  ts   <- grepl("^two_stage", new$row$method)
  pick <- function(mth, q) d[new$row$method == mth & new$row$quantity == q]
  dts  <- d[ts]
  check <- tibble::tibble(
    sim_id = new$meta$sim_id, scenario = new$meta$scenario,
    rel_maxd = max(d[new$row$method == "joint_4pl"]),
    abs_z_d = pick("joint_4pl_abs", "z"), abs_ct_d = pick("joint_4pl_abs", "CTmax_1hr"),
    two_stage_now_maxd = if (any(!is.na(dts))) max(dts, na.rm = TRUE) else NA_real_,
    kept_old = NROW(old$draws_abs), kept_new = NROW(new$draws_abs),
    uncrossed = new$meta$abs_uncrossed)
  # Rows: the re-run's joint rows and ALL of June's two-stage rows, in June's
  # order. Keying on the re-run's rows instead lost June rows the re-run did not
  # produce, and wrote all-NA rows for two-stage methods June did not have
  # (today's glmmTMB returns the t-interval beta-binomial on some datasets
  # where June's did not, and the reverse). Every row carries the re-run's
  # joint runtime, as run_one_sim() writes it.
  rows <- dplyr::bind_rows(new$row[!ts, , drop = FALSE],
                           old$row[grepl("^two_stage", old$row$method), , drop = FALSE])
  rows$runtime_sec <- new$row$runtime_sec[1]
  new$row <- rows
  new$meta$ts_bin_ok <- old$meta$ts_bin_ok
  new$meta$ts_bb_ok  <- old$meta$ts_bb_ok
  new$meta$two_stage <- "June 2026 run (kept)"
  # The fits reproduce exactly, but diagnose_tdt_fit() has since changed how it
  # counts treedepth hits (against the configured ceiling, not the run's own
  # maximum). Keep June's sampler diagnostics, so meta, diag and summary_conv
  # use one rule in every cell.
  dcols <- c("rhat_max", "ess_bulk_min", "ess_tail_min", "divergences", "treedepth_hits",
             "bfmi_min", "rhat_pass", "ess_pass", "divergence_pass", "treedepth_pass",
             "bfmi_pass", "all_pass")
  new$meta[dcols] <- old$meta[dcols]
  new$meta$treedepth_max <- NULL
  list(result = new, check = check)
}

#' Aggregate one scenario's raw files into the 8 objects the documents read
#'
#' @param label Scenario label.
#' @param raw_dir,out_dir Raw per-sim directory and output directory.
#' @return The scenario's MCSE summary (invisibly printed by the caller).
aggregate_scenario <- function(label, raw_dir, out_dir) {
  agg <- collect_raw(raw_dir)
  out <- function(prefix, obj)
    save_retry(obj, file.path(out_dir, sprintf("%s_%s.rds", prefix, label)))
  out("per_sim",      agg$per_sim)
  out("meta",         agg$meta)
  out("draws",        agg$draws)
  out("draws_abs",    agg$draws_abs)
  out("summary",      summarise_mcse(agg$per_sim))
  out("summary_conv", summarise_mcse(agg$per_sim, conv = TRUE, meta = agg$meta))
  out("diag",         diag_summary(agg$meta, label))
  out("diffs",        pairwise_diffs(agg$per_sim))
  summarise_mcse(agg$per_sim)
}

#' Bind all per-sim result files in a raw directory
#'
#' @param raw_dir Directory of `sim_#####.rds` files written by the runner.
#' @return List with `per_sim`, `meta`, `draws`, `draws_abs` (bound tibbles).
collect_raw <- function(raw_dir) {
  files <- list.files(raw_dir, pattern = "^sim_\\d+\\.rds$", full.names = TRUE)
  objs  <- lapply(files, readRDS)
  list(per_sim   = dplyr::bind_rows(lapply(objs, `[[`, "row")),
       meta      = dplyr::bind_rows(lapply(objs, `[[`, "meta")),
       draws     = dplyr::bind_rows(lapply(objs, `[[`, "draws")),
       draws_abs = dplyr::bind_rows(lapply(objs, `[[`, "draws_abs")))
}

#' Bias / coverage / width summary with Monte Carlo standard errors
#'
#' One row per (scenario, method, quantity). With `conv = TRUE`, restricts to
#' the simulations whose joint-4PL `all_pass` diagnostic is TRUE (same sims for
#' every method, so the comparison is on a common set) — the convergent-only
#' summary. Pass the run's `meta` table when `conv = TRUE`.
#'
#' @param per_sim Bound `per_sim` table from [collect_raw()].
#' @param conv Restrict to convergent joint-4PL sims? Default FALSE.
#' @param meta Meta table (required when `conv = TRUE`).
#' @return Tibble of summaries (rounded to 4 dp).
summarise_mcse <- function(per_sim, conv = FALSE, meta = NULL) {
  d <- dplyr::filter(per_sim, success)
  if (conv) {
    if (is.null(meta)) stop("summarise_mcse(conv = TRUE) needs `meta`.", call. = FALSE)
    # `joint_ok %in% TRUE` (not isTRUE(): that collapses the whole column to a
    # single value). Restrict to sims whose joint 4PL converged on all checks.
    ok <- meta |> dplyr::filter(joint_ok %in% TRUE, all_pass %in% TRUE) |>
      dplyr::pull(sim_id)
    d <- dplyr::filter(d, sim_id %in% ok)
  }
  # An empty `d` flows through to a correct 0-row summary (group_by + summarise
  # on no rows yields no groups), so no special-casing is needed.
  d |>
    dplyr::group_by(scenario, method, quantity) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean_bias = mean(bias, na.rm = TRUE),
      mcse_bias = stats::sd(bias, na.rm = TRUE) / sqrt(dplyr::n()),
      rmse = sqrt(mean(bias^2, na.rm = TRUE)),
      mcse_rmse = stats::sd(bias^2, na.rm = TRUE) /
                  (2 * sqrt(mean(bias^2, na.rm = TRUE)) * sqrt(dplyr::n())),
      coverage = mean(covered, na.rm = TRUE),
      mcse_cov = sqrt(mean(covered, na.rm = TRUE) *
                      (1 - mean(covered, na.rm = TRUE)) / dplyr::n()),
      med_width = stats::median(width, na.rm = TRUE),
      .groups = "drop") |>
    dplyr::mutate(dplyr::across(c(mean_bias, mcse_bias, rmse, mcse_rmse,
                                  coverage, mcse_cov, med_width), \(x) round(x, 4)))
}

#' Per-cell joint-4PL sampler diagnostic summary
#'
#' @param meta Meta table from [collect_raw()].
#' @param scenario Scenario label to stamp on the single summary row.
#' @return One-row tibble of median/worst-case diagnostics and pass rates.
diag_summary <- function(meta, scenario) {
  meta |>
    dplyr::filter(joint_ok) |>
    dplyr::summarise(
      n_total = dplyr::n(),
      n_joint_ok = sum(joint_ok, na.rm = TRUE),
      n_ts_bin_ok = sum(ts_bin_ok, na.rm = TRUE),
      n_ts_bb_ok = sum(ts_bb_ok, na.rm = TRUE),
      rhat_max_p50 = stats::median(rhat_max, na.rm = TRUE),
      rhat_max_p100 = max(rhat_max, na.rm = TRUE),
      ess_bulk_min_p50 = stats::median(ess_bulk_min, na.rm = TRUE),
      ess_bulk_min_p0 = suppressWarnings(min(ess_bulk_min, na.rm = TRUE)),
      ess_tail_min_p50 = stats::median(ess_tail_min, na.rm = TRUE),
      ess_tail_min_p0 = suppressWarnings(min(ess_tail_min, na.rm = TRUE)),
      bfmi_min_p50 = stats::median(bfmi_min, na.rm = TRUE),
      bfmi_min_p0 = suppressWarnings(min(bfmi_min, na.rm = TRUE)),
      n_div_total = sum(divergences, na.rm = TRUE),
      n_sims_with_div = sum(divergences > 0, na.rm = TRUE),
      n_treedepth_total = sum(treedepth_hits, na.rm = TRUE),
      n_sims_with_td = sum(treedepth_hits > 0, na.rm = TRUE),
      rhat_pass_rate = mean(rhat_pass, na.rm = TRUE),
      ess_pass_rate = mean(ess_pass, na.rm = TRUE),
      divergence_pass_rate = mean(divergence_pass, na.rm = TRUE),
      treedepth_pass_rate = mean(treedepth_pass, na.rm = TRUE),
      bfmi_pass_rate = mean(bfmi_pass, na.rm = TRUE),
      all_pass_rate = mean(all_pass, na.rm = TRUE),
      .groups = "drop") |>
    dplyr::mutate(scenario = scenario, .before = 1)
}

#' Pairwise within-simulation method differences
#'
#' All three contrasts (joint_4pl - two_stage_bin, joint_4pl - two_stage_bb,
#' two_stage_bb - two_stage_bin), aggregated across sims. Captures whether
#' methods agree on individual datasets — the within-sim contrast the
#' across-sim summary marginalises out.
#'
#' @param per_sim Bound `per_sim` table.
#' @return Tibble of per-contrast summaries (rounded to 4 dp).
pairwise_diffs <- function(per_sim) {
  wide <- per_sim |>
    dplyr::select(sim_id, scenario, method, quantity, estimate, success) |>
    tidyr::pivot_wider(names_from = method,
                       values_from = c(estimate, success), names_sep = "_")

  one <- function(m1, m2) {
    need <- c(paste0("estimate_", c(m1, m2)), paste0("success_", c(m1, m2)))
    if (!all(need %in% names(wide))) return(NULL)
    wide |>
      dplyr::filter(.data[[paste0("success_", m1)]],
                    .data[[paste0("success_", m2)]]) |>
      dplyr::transmute(sim_id, scenario, quantity, m1 = m1, m2 = m2,
                       diff = .data[[paste0("estimate_", m1)]] -
                              .data[[paste0("estimate_", m2)]])
  }
  diffs <- dplyr::bind_rows(one("joint_4pl", "two_stage_bin"),
                            one("joint_4pl", "two_stage_bb"),
                            one("two_stage_bb", "two_stage_bin"))
  if (nrow(diffs) == 0L) return(diffs)

  diffs |>
    dplyr::group_by(scenario, m1, m2, quantity) |>
    dplyr::summarise(n = dplyr::n(),
      mean_diff = mean(diff, na.rm = TRUE),
      mcse_diff = stats::sd(diff, na.rm = TRUE) / sqrt(dplyr::n()),
      median_diff = stats::median(diff, na.rm = TRUE),
      diff_q025 = stats::quantile(diff, 0.025, na.rm = TRUE, names = FALSE),
      diff_q975 = stats::quantile(diff, 0.975, na.rm = TRUE, names = FALSE),
      .groups = "drop") |>
    dplyr::mutate(dplyr::across(c(mean_diff, mcse_diff, median_diff,
                                  diff_q025, diff_q975), \(x) round(x, 4)))
}

# ----------------------------------------------------------------------------
# 7. Small utilities used by the runner
# ----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

#' NA -> NULL (so a tribble's NA cell becomes "keep the DGP preset")
na_to_null <- function(x) if (length(x) == 1 && is.na(x)) NULL else x

#' Simulation ids already complete in a raw directory
#'
#' @param raw_dir Directory of `sim_#####.rds` files.
#' @param force If TRUE, treat nothing as done (forces a full refit).
#' @return Integer vector of completed sim ids.
done_ids <- function(raw_dir, force = FALSE) {
  if (force) return(integer(0))
  f <- list.files(raw_dir, pattern = "^sim_\\d+\\.rds$")
  if (!length(f)) integer(0) else as.integer(sub("^sim_(\\d+)\\.rds$", "\\1", f))
}

#' Run a function over ids on a PSOCK cluster
#'
#' PSOCK (not `mclapply`'s fork) because cmdstanr launches Stan via processx,
#' which is not fork-safe on macOS. Workers load bayesTLS and source this file,
#' then `fun` runs with the exported objects visible.
#'
#' Robust saveRDS: retry on transient I/O errors.
#'
#' Dropbox/CloudStorage and other networked file providers occasionally fail a
#' write with "cannot open the connection" / "Interrupted system call" (EINTR)
#' when the provider is busy. A single such failure crashed an entire run
#' (2026-06-13) at the final aggregate save. This retries the write a few times
#' with a pause so a transient hiccup does not lose hours of compute.
#'
#' @param object,file Passed to [saveRDS()].
#' @param tries,wait Max attempts and seconds to wait between them.
save_retry <- function(object, file, tries = 10, wait = 6) {
  for (i in seq_len(tries)) {
    ok <- tryCatch({ saveRDS(object, file); TRUE },
                   error = function(e) {
                     message(sprintf("  saveRDS retry %d/%d for %s: %s",
                                     i, tries, basename(file), conditionMessage(e)))
                     FALSE
                   })
    if (ok) return(invisible(file))
    Sys.sleep(wait)
  }
  stop("save_retry: gave up after ", tries, " attempts writing ", file)
}

#' @param ids Integer ids to map over.
#' @param fun Function of one id (the per-sim worker).
#' @param workers Cluster size; 1 runs serially (easier debugging).
#' @param exports Names of objects in `envir` the workers need.
#' @param envir Environment to export from. Default the caller's.
#' @param lib Optional library path put first on each worker's `.libPaths()`.
run_parallel <- function(ids, fun, workers, exports, envir = parent.frame(),
                         lib = NULL,
                         src = here::here("scripts", "simulations", "sim_functions.R")) {
  if (workers <= 1L) { invisible(lapply(ids, fun)); return(invisible()) }
  cl <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  # An optional private library (e.g. a frozen bayesTLS) goes first on each
  # worker's library path, as in the master, before bayesTLS is attached.
  parallel::clusterCall(cl, function(l) if (!is.null(l)) .libPaths(c(l, .libPaths())), lib)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages({
    library(bayesTLS); library(dplyr); library(tibble); library(tidyr)
  }))
  # Workers source `src` (a long run can pass a frozen copy, so that editing
  # this file mid-run cannot change the cells still to come).
  parallel::clusterCall(cl, function(f) source(f), normalizePath(src))
  parallel::clusterExport(cl, exports, envir = envir)
  invisible(parallel::parLapply(cl, ids, fun))
}

#' Pre-flight every scenario before any long compute
#'
#' For each scenario: build its truth (the DGP feasibility guard in
#' [compute_ols_truth()] fires on infeasible u/ell/k), generate one dataset on a
#' fixed seed (catches NaN survival), and hash the survival vector. Any two
#' scenarios with different arguments must hash differently — a duplicate means
#' a flag isn't threaded through to data generation (the 2026-05-15 bug class).
#'
#' @param scenarios Data frame of scenario rows (label + DGP args).
#' @return Invisibly, the hash table; stops on infeasible DGP or duplicate hash.
preflight <- function(scenarios) {
  # A scenario's data is fully determined by its *resolved* generating spec:
  # the family, replication, design grid, and the 4PL truth parameters AFTER
  # presets + overrides are applied — not the raw scenario arguments. Two
  # scenarios that resolve to the same spec legitimately share a preflight hash
  # (e.g. "full"/"tmax405" expand to the same grid; the asym_u DGP equals a
  # baseline DGP with u_beta1 = -0.01). Two with *different* resolved specs must
  # not — that would mean a flag never reached data generation.
  spec_key <- function(truth, n_reps) {
    sp <- design_spec(truth$design)
    paste(truth$family, n_reps,
          paste(sp$temps, collapse = ","), paste(sp$durations, collapse = ","),
          truth$u, truth$ell, truth$k, truth$m_beta0, truth$m_beta1,
          truth$u_beta1, truth$ell_beta1, truth$k_beta1, truth$T_bar, truth$phi,
          sep = "|")
  }

  hashes <- character(nrow(scenarios))
  keys   <- character(nrow(scenarios))
  for (i in seq_len(nrow(scenarios))) {
    sc <- scenarios[i, ]
    truth <- sim_truth(dgp = sc$dgp, family = sc$family, design = sc$design,
                       u_0 = na_to_null(sc$u_0), ell_0 = na_to_null(sc$ell_0),
                       u_beta1 = na_to_null(sc$u_beta1))
    d <- sim_dataset(n_reps = sc$n_reps, seed = 1L, truth = truth)
    if (any(is.na(d$y)))
      stop("preflight: ", sc$label, " produced NA survival — DGP infeasible.",
           call. = FALSE)
    hashes[i] <- substr(rlang::hash(d$y), 1, 12)
    keys[i]   <- spec_key(truth, sc$n_reps)
    message(sprintf("  OK %-22s z_true=%.3f  CTmax_true=%.3f  hash=%s",
                    sc$label, truth$z_true, truth$CTmax_1hr_true, hashes[i]))
  }

  dup <- hashes[duplicated(hashes)]
  if (length(dup)) {
    for (h in unique(dup)) {
      labs <- scenarios$label[hashes == h]
      if (length(unique(keys[hashes == h])) > 1L)
        stop("preflight: scenarios with different specs share data hash ",
             h, ": ", paste(labs, collapse = ", "),
             " — a flag isn't reaching data generation.", call. = FALSE)
      message(sprintf("  NOTE expected duplicate hash %s: %s (identical spec).",
                      h, paste(labs, collapse = ", ")))
    }
  }
  invisible(hashes)
}
