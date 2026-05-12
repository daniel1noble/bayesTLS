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
