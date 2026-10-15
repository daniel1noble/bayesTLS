# Daniel's proposed formulation (2026-09-26), tested against the same DGP as
# notes/2026-09-26-naive-vs-conditional-sim.R (same seeds 101:120, same truth):
#
#   rows  = individual x check day, y = 0/1 CUMULATIVE dead status (dead stays
#           dead; post-death rows included, as literally proposed)
#   mean  = 1 - S(t) with the 4PL itself (marginal, NOT conditional)
#   D1: no random effect            [prediction: identical posterior to the
#                                    cell-level naive binomial, by aggregation]
#   D2: + (1 | id) on mid           [the proposal: does the RE repair the
#                                    absorbing-state dependence?]
#   D3: + (1 + lt | id) on mid      [random time slope; 3 reps only — expect
#                                    weak identification from one event/animal]
#
# Mid-parameterisation used for comparability with the other arms; the package
# audit (2026-09-25) verified direct CTmax/z emits the same model (z = -1/b1).

suppressPackageStartupMessages({ library(brms) })
options(mc.cores = 1)

S <- Sys.getenv("SCRATCH", unset = dirname(sub("--file=", "",
       grep("--file=", commandArgs(FALSE), value = TRUE)[1])))

temps      <- c(37, 39, 41, 43)
t_mean     <- mean(temps)
n_per_temp <- 40
check_days <- 14
true_b0    <- 0.5
true_z     <- 4
true_b1    <- -1 / true_z
true_k     <- 8
seeds      <- 101:120

sim_cohort <- function(seed) {
  set.seed(seed)
  do.call(rbind, lapply(seq_along(temps), function(i) {
    temp_c <- temps[i] - t_mean
    mid    <- true_b0 + true_b1 * temp_c
    U      <- runif(n_per_temp)
    data.frame(temp_c = temp_c,
               t_dead = 10^(mid + log(1 / U - 1) / true_k))
  }))
}

## Daniel's rows: individual x every check day, cumulative dead status.
make_indiv_cum <- function(sim) {
  sim$id <- seq_len(nrow(sim))
  out <- do.call(rbind, lapply(seq_len(nrow(sim)), function(r) {
    j <- seq_len(check_days)
    data.frame(id = sim$id[r], temp_c = sim$temp_c[r],
               lt = log10(j), dead = as.integer(sim$t_dead[r] <= j))
  }))
  stopifnot(all(tapply(out$dead, out$id, function(x) all(diff(x) >= 0))))
  out
}

form_d1 <- bf(dead ~ 1 - 1 / (1 + exp(exp(logk) * (lt - mid))),
              logk ~ 1, mid ~ temp_c, nl = TRUE)
form_d2 <- bf(dead ~ 1 - 1 / (1 + exp(exp(logk) * (lt - mid))),
              logk ~ 1, mid ~ temp_c + (1 | id), nl = TRUE)
form_d3 <- bf(dead ~ 1 - 1 / (1 + exp(exp(logk) * (lt - mid))),
              logk ~ 1, mid ~ temp_c + (1 + lt | id), nl = TRUE)

priors_mean <- c(prior(normal(1.5, 1), nlpar = "logk"),
                 prior(normal(0.5, 1), nlpar = "mid", coef = "Intercept"),
                 prior(normal(0, 0.5), nlpar = "mid", coef = "temp_c"))
prior_sd  <- prior(exponential(2), class = "sd", nlpar = "mid")
prior_cor <- prior(lkj(2), class = "cor")
init_fun  <- function() list(b_logk = as.array(1.5), b_mid = c(0.5, 0))

specs <- list(
  D1 = list(formula = form_d1, prior = priors_mean,                          reps = 1:20),
  D2 = list(formula = form_d2, prior = c(priors_mean, prior_sd),             reps = 1:20),
  D3 = list(formula = form_d3, prior = c(priors_mean, prior_sd, prior_cor),  reps = 1:3))

fit_one <- function(base, sp, dat, seed) {
  if (is.null(base))
    brm(sp$formula, data = dat, family = bernoulli(link = "identity"),
        prior = sp$prior, chains = 1, iter = 800, warmup = 400,
        init = init_fun, seed = seed, refresh = 0, backend = "cmdstanr",
        silent = 2)
  else
    update(base, newdata = dat, recompile = FALSE, chains = 1, iter = 800,
           warmup = 400, init = init_fun, seed = seed, refresh = 0, silent = 2)
}

res <- list(); base_fits <- list()
for (m in names(specs)) {
  sp <- specs[[m]]
  for (r in sp$reps) {
    seed <- seeds[r]
    dat  <- make_indiv_cum(sim_cohort(seed))
    row <- data.frame(rep = r, seed = seed, method = m, fit_ok = 0L,
                      z_med = NA, z_lo = NA, z_hi = NA, width = NA,
                      covers = NA, n_div = NA, max_rhat = NA,
                      sd_mid_med = NA, sd_lt_med = NA, err = "")
    t0 <- Sys.time()
    fit <- tryCatch(fit_one(base_fits[[m]], sp, dat, seed),
                    error = function(e) e)
    if (inherits(fit, "brmsfit")) {
      if (is.null(base_fits[[m]])) base_fits[[m]] <- fit
      dr  <- as.data.frame(fit, variable = "b_mid_temp_c")
      b1  <- dr[[1]]
      z   <- -1 / b1
      qs  <- unname(quantile(z[b1 < 0], c(.025, .5, .975)))
      sm  <- summary(fit)
      row$fit_ok <- 1L
      row$z_med <- qs[2]; row$z_lo <- qs[1]; row$z_hi <- qs[3]
      row$width <- qs[3] - qs[1]
      row$covers <- as.integer(qs[1] <= true_z && true_z <= qs[3])
      row$n_div <- sum(brms::nuts_params(fit, pars = "divergent__")$Value)
      row$max_rhat <- max(brms::rhat(fit), na.rm = TRUE)
      sds <- tryCatch(as.data.frame(fit, variable = "sd_id__mid_Intercept"),
                      error = function(e) NULL)
      if (!is.null(sds)) row$sd_mid_med <- median(sds[[1]])
      sdl <- tryCatch(as.data.frame(fit, variable = "sd_id__mid_lt"),
                      error = function(e) NULL)
      if (!is.null(sdl)) row$sd_lt_med <- median(sdl[[1]])
    } else row$err <- conditionMessage(fit)
    cat(sprintf("[%s rep %02d] ok=%d z=%.3f width=%.3f covers=%s div=%s rhat=%.3f sd=%.3f (%.1fs)\n",
                m, r, row$fit_ok, row$z_med, row$width, row$covers, row$n_div,
                row$max_rhat, row$sd_mid_med,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    res[[length(res) + 1]] <- row
  }
}
res <- do.call(rbind, res)
write.csv(res, file.path(S, "daniel_formulation_results.csv"), row.names = FALSE)

cat("\n--- summary (truth z = 4) ---\n")
for (m in names(specs)) {
  x <- res[res$method == m & res$fit_ok == 1, ]
  cat(sprintf("%s  n=%d  z_med(mean)=%.3f  width(mean)=%.3f  coverage=%d/%d  div_total=%d  max_rhat=%.3f  sd_mid(med)=%.3f  sd_lt(med)=%.3f\n",
              m, nrow(x), mean(x$z_med), mean(x$width), sum(x$covers), nrow(x),
              sum(x$n_div), max(x$max_rhat), median(x$sd_mid_med),
              median(x$sd_lt_med)))
}
cat("\nReference (same seeds, from naive_vs_conditional_results.csv):\n")
prev <- read.csv(file.path(S, "naive_vs_conditional_results.csv"))
for (m in c("A", "B", "C")) {
  x <- prev[prev$method == m, ]
  cat(sprintf("%s  width(mean)=%.3f  coverage=%d/20\n", m, mean(x$width), sum(x$covers)))
}
