# Reusable simulator for tests. Produces a beta-binomial-noised survival
# dataset with known truth, ready to feed to standardize_data().

simulate_tdt <- function(temps           = c(30, 33, 36),
                         durations_hours = c(0.05, 0.5, 5, 50),
                         n_rep           = 5,
                         n_per_rep       = 30,
                         phi             = 8,
                         truth           = list(ell     = 0.05,
                                                u       = 0.95,
                                                k       = 6,
                                                m_beta0 = 1.5,
                                                m_beta1 = -0.18,
                                                T_bar   = 33),
                         seed            = 1) {
  set.seed(seed)
  d <- expand.grid(T = temps, t_hours = durations_hours, rep = seq_len(n_rep))
  d$log10_t_min <- log10(d$t_hours * 60)
  d$mid_t       <- truth$m_beta0 + truth$m_beta1 * (d$T - truth$T_bar)
  d$p_true      <- truth$ell + (truth$u - truth$ell) /
                   (1 + exp(truth$k * (d$log10_t_min - d$mid_t)))
  d$n_trials    <- n_per_rep
  d$p_draw      <- rbeta(nrow(d), d$p_true * phi, (1 - d$p_true) * phi)
  d$y_alive     <- rbinom(nrow(d), size = d$n_trials, prob = d$p_draw)
  attr(d, "truth") <- truth
  d
}

# Continuous-proportion variant: draws the response directly from a Beta whose
# mean is the 4PL p_true (no binomial layer), for testing the Beta family path.
# Same truth/design as simulate_tdt() so truth_summary() applies unchanged.
simulate_tdt_beta <- function(temps           = c(30, 33, 36),
                              durations_hours = c(0.05, 0.5, 5, 50),
                              n_rep           = 8,
                              phi             = 8,
                              truth           = list(ell     = 0.05,
                                                     u       = 0.95,
                                                     k       = 6,
                                                     m_beta0 = 1.5,
                                                     m_beta1 = -0.18,
                                                     T_bar   = 33),
                              n_zeros         = 0,
                              seed            = 1) {
  set.seed(seed)
  d <- expand.grid(T = temps, t_hours = durations_hours, rep = seq_len(n_rep))
  d$log10_t_min <- log10(d$t_hours * 60)
  d$mid_t       <- truth$m_beta0 + truth$m_beta1 * (d$T - truth$T_bar)
  d$p_true      <- truth$ell + (truth$u - truth$ell) /
                   (1 + exp(truth$k * (d$log10_t_min - d$mid_t)))
  d$p_true      <- pmin(pmax(d$p_true, 1e-6), 1 - 1e-6)
  d$y_prop      <- rbeta(nrow(d), d$p_true * phi, (1 - d$p_true) * phi)
  if (n_zeros > 0) {
    hot <- which(d$T == max(temps) &
                 d$t_hours >= stats::median(durations_hours))
    if (length(hot) > 0) {
      d$y_prop[sample(hot, min(n_zeros, length(hot)))] <- 0
    }
  }
  attr(d, "truth") <- truth
  d
}

# Truth values derived from the simulator's default parameters.
truth_summary <- function(truth = list(ell = 0.05, u = 0.95, k = 6,
                                       m_beta0 = 1.5, m_beta1 = -0.18,
                                       T_bar = 33)) {
  z       <- -1 / truth$m_beta1
  alpha   <- truth$m_beta0 +
             (1 / truth$k) * log((truth$u - 0.5) / (0.5 - truth$ell)) -
             truth$m_beta1 * truth$T_bar
  CT1     <- (log10(60) - alpha) / truth$m_beta1
  list(z = z, CTmax_1hr = CT1, truth = truth)
}

# Simulator with KNOWN temperature SLOPES on the asymptotes/steepness, not just
# the midpoint. Used to test that extract_4pl_pars() faithfully recovers the
# per-temperature low/up/k/mid (the slopes the old constant-shape extractor
# discarded). The natural-scale sub-parameters are linear in temp_c (centred at
# T_bar): low ~ 1 (constant), up ~ temp_c (declines), k ~ temp_c (steepens),
# mid ~ temp_c (declines). logd is on the log10(hours) scale to match
# standardize_data(duration_unit = "hours") (which sets logd = log10(duration)).
simulate_tdt_shape <- function(temps    = c(30, 31.5, 33, 34.5, 36),
                               durs_h   = c(0.05, 0.2, 0.8, 3, 12, 48),
                               n_rep    = 8,
                               n_per    = 40,
                               phi      = 12,
                               T_bar    = 33,
                               truth    = list(low0 = 0.05,
                                               up0 = 0.90, up_slope = -0.02,
                                               logk0 = log(6), logk_slope = 0.10,
                                               mid0 = 1.4, mid_slope = -0.18),
                               seed     = 101) {
  set.seed(seed)
  d  <- expand.grid(T = temps, t_hours = durs_h, rep = seq_len(n_rep))
  tc <- d$T - T_bar
  logd <- log10(d$t_hours)
  low <- rep(truth$low0, nrow(d))
  up  <- truth$up0 + truth$up_slope * tc
  k   <- exp(truth$logk0 + truth$logk_slope * tc)
  mid <- truth$mid0 + truth$mid_slope * tc
  d$p_true   <- pmin(pmax(low + (up - low) / (1 + exp(k * (logd - mid))),
                          1e-5), 1 - 1e-5)
  d$n_trials <- n_per
  d$p_draw   <- rbeta(nrow(d), d$p_true * phi, (1 - d$p_true) * phi)
  d$y_alive  <- rbinom(nrow(d), d$n_trials, d$p_draw)
  attr(d, "truth_shape") <- c(truth, list(T_bar = T_bar, temps = temps))
  d
}

# Natural-scale truth for simulate_tdt_shape() at a vector of (uncentred) temps.
truth_shape_at <- function(temps, ts) {
  tc <- temps - ts$T_bar
  data.frame(temp = temps,
             low  = rep(ts$low0, length(temps)),
             up   = ts$up0 + ts$up_slope * tc,
             k    = exp(ts$logk0 + ts$logk_slope * tc),
             mid  = ts$mid0 + ts$mid_slope * tc)
}

# Two-group simulator with DELIBERATELY different per-group z (large gap so a
# swapped group mapping can't pass by coincidence): group A m_beta1 = -0.18
# (z = 5.56), group B m_beta1 = -0.30 (z = 3.33). Same design otherwise. Used for
# the grouped/per-group gated tests (fit with ctmax/z ~ 0 + grp).
simulate_tdt_2group <- function(seed = 20260622) {
  tA <- list(ell = 0.05, u = 0.95, k = 6, m_beta0 = 1.5, m_beta1 = -0.18, T_bar = 33)
  tB <- list(ell = 0.05, u = 0.95, k = 6, m_beta0 = 1.2, m_beta1 = -0.30, T_bar = 33)
  a <- simulate_tdt(seed = seed,     truth = tA); a$grp <- "A"
  b <- simulate_tdt(seed = seed + 1, truth = tB); b$grp <- "B"
  d <- rbind(a, b); d$grp <- factor(d$grp)
  attr(d, "truth_2group") <- list(A = truth_summary(tA), B = truth_summary(tB))
  d
}

# PURE-BINOMIAL two-group simulator (no beta-binomial overdispersion: y_alive is
# drawn straight from Binomial(n, p_true)). Distinct per-group z AND CTmax (gap
# big enough that a label swap can't pass), with a small per-group asymptote
# difference (group B has a lower upper asymptote) so per-group `up` recovery is
# also exercised. Fit with family = binomial(link = "identity").
simulate_tdt_2group_binom <- function(temps = c(29, 32, 35, 38),
                                      durations_hours = c(0.05, 0.5, 5, 50),
                                      n_rep = 12, n_per_rep = 45, seed = 20260626) {
  tA <- list(ell = 0.03, u = 0.97, k = 6, m_beta0 = 1.5, m_beta1 = -0.16, T_bar = 33)
  tB <- list(ell = 0.03, u = 0.90, k = 6, m_beta0 = 1.1, m_beta1 = -0.28, T_bar = 33)
  one <- function(tr, sd, grp) {
    set.seed(sd)
    d <- expand.grid(T = temps, t_hours = durations_hours, rep = seq_len(n_rep))
    d$log10_t_min <- log10(d$t_hours * 60)
    d$mid_t  <- tr$m_beta0 + tr$m_beta1 * (d$T - tr$T_bar)
    d$p_true <- tr$ell + (tr$u - tr$ell) / (1 + exp(tr$k * (d$log10_t_min - d$mid_t)))
    d$n_trials <- n_per_rep
    d$y_alive  <- rbinom(nrow(d), d$n_trials, d$p_true)   # PURE binomial, no beta layer
    d$grp <- grp
    d
  }
  d <- rbind(one(tA, seed, "A"), one(tB, seed + 1, "B"))
  d$grp <- factor(d$grp)
  attr(d, "truth_2group") <- list(A = c(truth_summary(tA), list(low = tA$ell, up = tA$u)),
                                  B = c(truth_summary(tB), list(low = tB$ell, up = tB$u)))
  d
}

# Continuous-proportion (Beta) two-group simulator: response drawn directly from
# a Beta whose mean is the 4PL p_true (no binomial layer), with distinct
# per-group z and CTmax. Fit with standardize_data(proportion=) ->
# fit_4pl(family = Beta(link = "identity")).
simulate_tdt_beta_2group <- function(temps = c(30, 33, 36),
                                     durations_hours = c(0.05, 0.5, 5, 50),
                                     n_rep = 12, phi = 12, seed = 20260626) {
  tA <- list(ell = 0.03, u = 0.97, k = 6, m_beta0 = 1.5, m_beta1 = -0.16, T_bar = 33)
  tB <- list(ell = 0.03, u = 0.97, k = 6, m_beta0 = 1.1, m_beta1 = -0.28, T_bar = 33)
  one <- function(tr, sd, grp) {
    set.seed(sd)
    d <- expand.grid(T = temps, t_hours = durations_hours, rep = seq_len(n_rep))
    d$log10_t_min <- log10(d$t_hours * 60)
    d$mid_t  <- tr$m_beta0 + tr$m_beta1 * (d$T - tr$T_bar)
    d$p_true <- tr$ell + (tr$u - tr$ell) / (1 + exp(tr$k * (d$log10_t_min - d$mid_t)))
    d$p_true <- pmin(pmax(d$p_true, 1e-6), 1 - 1e-6)
    d$y_prop <- rbeta(nrow(d), d$p_true * phi, (1 - d$p_true) * phi)
    d$grp <- grp
    d
  }
  d <- rbind(one(tA, seed, "A"), one(tB, seed + 1, "B"))
  d$grp <- factor(d$grp)
  attr(d, "truth_2group") <- list(A = truth_summary(tA), B = truth_summary(tB))
  d
}

# Two interacting CONTINUOUS moderators (x1, x2), CENTRED at 0. Plants the direct
# parameterisation surfaces:
#   CTmaxdev(x1,x2) = a0 + a1 x1 + a2 x2 + a3 x1:x2   (deg C, vs T_bar)
#   logz(x1,x2)     = z0 + z1 x1 + z2 x2 + z3 x1:x2   (z = exp(logz))
# The relative-direct backbone the model fits is
#   mid(T) = log10_tref - (temp_c - CTmaxdev)/exp(logz)
# with log10_tref = 0 on the HOURS axis (t_ref = 60 min = 1 h, logd = log10(hours)).
# `a3`/`z3` are the x1:x2 interaction coefficients; set them to 0 for the null
# companion fit used in the interaction-detectability check. Asymptotes constant.
simulate_tdt_cont2 <- function(a   = c(a0 = 1.0, a1 = 0.8, a2 = -0.5, a3 = 0.6),
                               zc  = c(z0 = log(5), z1 = 0.15, z2 = -0.10, z3 = 0.20),
                               T_bar = 35,
                               temps  = c(32, 34, 36, 38),
                               durs_h = c(0.1, 0.5, 2, 10),
                               x_levels = c(-1, -0.5, 0, 0.5, 1),
                               n_rep = 6, n_per = 30,
                               low = 0.03, up = 0.97, k = 6,
                               seed = 303) {
  set.seed(seed)
  covpts <- expand.grid(x1 = x_levels, x2 = x_levels)
  d <- do.call(rbind, lapply(seq_len(nrow(covpts)), function(i) {
    g <- expand.grid(T = temps, t_hours = durs_h)
    g$x1 <- covpts$x1[i]; g$x2 <- covpts$x2[i]; g
  }))
  d <- d[rep(seq_len(nrow(d)), each = n_rep), ]
  tc    <- d$T - T_bar
  CTdev <- a["a0"] + a["a1"] * d$x1 + a["a2"] * d$x2 + a["a3"] * d$x1 * d$x2
  logz  <- zc["z0"] + zc["z1"] * d$x1 + zc["z2"] * d$x2 + zc["z3"] * d$x1 * d$x2
  z     <- exp(logz)
  mid   <- 0 - (tc - CTdev) / z              # log10_tref(1h) = 0 on the hours axis
  # Local log10(hours) for the DGP only; do NOT name it `logd` (standardize_data
  # reserves that name and would warn about overwriting it).
  log10_t_h <- log10(d$t_hours)
  d$p_true <- pmin(pmax(low + (up - low) / (1 + exp(k * (log10_t_h - mid))), 1e-5), 1 - 1e-5)
  d$n_trials <- n_per
  d$y_alive  <- rbinom(nrow(d), d$n_trials, d$p_true)
  attr(d, "truth_cont2") <- list(a = a, zc = zc, T_bar = T_bar,
                                 low = low, up = up, k = k)
  d
}
