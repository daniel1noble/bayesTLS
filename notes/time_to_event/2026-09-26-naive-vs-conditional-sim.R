#!/usr/bin/env Rscript
# ==============================================================================
# naive_vs_conditional_sim.R
#
# QUESTION: a collaborator monitors the SAME cohort at daily checkpoints and
# proposes feeding CUMULATIVE survival proportions at each checkpoint into the
# bayesTLS batch-binomial model as if they were independent cells. How
# overconfident is that naive reformat in practice, and does beta-binomial
# overdispersion (the package default family) partially rescue it?
#
# Truth: individual log-logistic time-to-death implied by the bayesTLS 4PL
# (low = 0, up = 1): S(t; T) = 1 / (1 + exp(k (log10 t - mid))),
# mid = b0 + b1 * temp_c, z = -1/b1. (Mappings verified in tte_math_check.R.)
#
# Three fits per replicate, all with the SAME 4PL mean structure and the SAME
# weakly-informative priors on logk / mid (as in tte_prototype.R):
#   (A) NAIVE-BINOMIAL      rows = temp x checkpoint (4 x 14 = 56),
#                           alive_j | trials(40), mean S(t_j)  [cumulative --
#                           exactly what the current package would be handed]
#   (B) NAIVE-BETABINOMIAL  same 56 rows, beta_binomial(identity), phi
#                           estimated with the package default prior
#                           gamma(2, 0.1)  [the package default family]
#   (C) CORRECT             person-period conditional Bernoulli,
#                           p_j = 1 - S(t_j)/S(t_{j-1})  [tte_prototype.R
#                           construction; log1p_exp spelling of the SAME
#                           ratio -- algebraically identical, verified in
#                           tte_math_check.R claim (b), numerically stable
#                           across seeds]
#
# Design per replicate: temps 37/39/41/43 C (centred at 40), 40 individuals
# per temp, true log10 t50 = 0.5 - 0.25 * temp_c (z = 4 C), k = 8, daily
# checks days 1..14, right-censored at day 14. R = 20 replicates, seeds
# 101..120 (seed drives BOTH data generation and the sampler).
#
# Sampling: cmdstanr backend, 1 chain, iter 800 (400 warmup), refresh = 0,
# NO file caching. Each model is compiled ONCE (rep 1) and reused via
# update(fit, newdata = ...) for reps 2..20.
#
# Outputs:
#   naive_vs_conditional_results.csv  one row per replicate x method
#   naive_vs_conditional_summary.txt  aggregated report
# ==============================================================================

SCRATCH <- "/private/tmp/claude-501/-Users-noble-Library-CloudStorage-Dropbox-1-Research-1-Manuscripts-In-Preparation-tls-model-equivalence/5004b512-34de-4b37-8d53-7208be821650/scratchpad"
CSV_OUT <- file.path(SCRATCH, "naive_vs_conditional_results.csv")
TXT_OUT <- file.path(SCRATCH, "naive_vs_conditional_summary.txt")

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
})
options(brms.backend = "cmdstanr")

## ---- fixed design ------------------------------------------------------------
temps      <- c(37, 39, 41, 43)
n_per_temp <- 40
t_mean     <- mean(temps)            # centring constant (40 C)
check_days <- 14
true_b0    <- 0.5
true_z     <- 4
true_b1    <- -1 / true_z            # -0.25
true_k     <- 8
seeds      <- 101:120
R          <- length(seeds)

## ---- data generation ---------------------------------------------------------
## Individual death times from the log-logistic implied by the 4PL:
## S(t) = 1/(1+exp(k(log10 t - mid))) = U  =>  log10 t = mid + log(1/U - 1)/k
sim_cohort <- function(seed) {
  set.seed(seed)
  do.call(rbind, lapply(seq_along(temps), function(i) {
    temp_c <- temps[i] - t_mean
    mid    <- true_b0 + true_b1 * temp_c
    U      <- runif(n_per_temp)
    data.frame(
      temp   = temps[i],
      temp_c = temp_c,
      t_dead = 10^(mid + log(1 / U - 1) / true_k)
    )
  }))
}

## NAIVE reformat: per-treatment cumulative survival at every checkpoint,
## one row per temp x checkpoint, all with trials = the starting n (40).
## Alive at check j  <=>  t_dead > j (death in (j-1, j] counted dead at j).
make_naive <- function(sim) {
  out <- do.call(rbind, lapply(temps, function(Tt) {
    td <- sim$t_dead[sim$temp == Tt]
    data.frame(
      temp_c = Tt - t_mean,
      day    = seq_len(check_days),
      lt     = log10(seq_len(check_days)),
      alive  = vapply(seq_len(check_days),
                      function(j) sum(td > j), integer(1)),
      n0     = n_per_temp
    )
  }))
  ## sanity: cumulative survival is non-increasing within each temperature
  for (Tt in temps) {
    a <- out$alive[out$temp_c == Tt - t_mean]
    stopifnot("naive alive counts non-increasing in day" = all(diff(a) <= 0))
  }
  out
}

## CORRECT person-period restructure: one row per individual x interval while
## alive at interval start; death indicator for that interval.
make_pp <- function(sim) {
  death_int <- ifelse(sim$t_dead <= check_days, ceiling(sim$t_dead), NA)
  do.call(rbind, lapply(seq_len(nrow(sim)), function(r) {
    n_int <- if (is.na(death_int[r])) check_days else death_int[r]
    j     <- seq_len(n_int)
    data.frame(
      temp_c = sim$temp_c[r],
      ## -10 stands in for log10(0) = -Inf: S(10^-10 d) = 1 to machine
      ## precision for any plausible (k, mid); keeps the design finite.
      lt_lo  = ifelse(j == 1, -10, log10(j - 1)),
      lt_hi  = log10(j),
      death  = as.integer(!is.na(death_int[r]) & j == n_int)
    )
  }))
}

## ---- formulas, priors, inits -------------------------------------------------
form_naive <- bf(
  alive | trials(n0) ~ 1 / (1 + exp(exp(logk) * (lt - mid))),
  logk ~ 1,
  mid  ~ temp_c,
  nl = TRUE
)
form_pp <- bf(
  death ~ 1 - exp(log1p_exp(exp(logk) * (lt_lo - mid)) -
                  log1p_exp(exp(logk) * (lt_hi - mid))),
  logk ~ 1,
  mid  ~ temp_c,
  nl = TRUE
)

priors_mean <- c(
  prior(normal(1.5, 1), nlpar = "logk"),
  prior(normal(0.5, 1), nlpar = "mid", coef = "Intercept"),
  prior(normal(0, 0.5), nlpar = "mid", coef = "temp_c")
)
## Package default phi prior (R/priors.R): gamma(2, 0.1)
priors_bb <- c(priors_mean, prior(gamma(2, 0.1), class = "phi"))

## Prior-centre inits (identity-link likelihoods reject regions where the
## mean underflows, so start each chain somewhere sane -- as in the prototype).
init_mean <- function() list(b_logk = as.array(1.5), b_mid = c(0.5, 0))
init_bb   <- function() list(b_logk = as.array(1.5), b_mid = c(0.5, 0),
                             phi = 20)

## ---- metric extraction -------------------------------------------------------
extract_metrics <- function(fit, has_phi) {
  draws <- posterior::as_draws_df(fit)
  b1    <- draws$b_mid_temp_c
  z     <- -1 / b1
  zq    <- unname(quantile(z, c(0.025, 0.5, 0.975)))
  vars  <- c("b_logk_Intercept", "b_mid_Intercept", "b_mid_temp_c",
             if (has_phi) "phi")
  st    <- posterior::summarise_draws(
    posterior::subset_draws(draws, variable = vars)
  )
  np    <- brms::nuts_params(fit)
  list(
    z_med        = zq[2],
    z_lo         = zq[1],
    z_hi         = zq[3],
    width        = zq[3] - zq[1],
    covers       = as.integer(zq[1] <= true_z && true_z <= zq[3]),
    n_div        = sum(np$Value[np$Parameter == "divergent__"]),
    max_rhat     = max(st$rhat, na.rm = TRUE),
    min_ess      = min(st$ess_bulk, na.rm = TRUE),
    n_b1_nonneg  = sum(b1 >= 0),     # draws where z = -1/b1 is undefined-signed
    phi_med      = if (has_phi) median(draws$phi) else NA_real_
  )
}

row_template <- function(rep, seed, method) {
  data.frame(rep = rep, seed = seed, method = method, fit_ok = 0L,
             z_med = NA_real_, z_lo = NA_real_, z_hi = NA_real_,
             width = NA_real_, covers = NA_integer_, n_div = NA_real_,
             max_rhat = NA_real_, min_ess = NA_real_,
             n_b1_nonneg = NA_integer_, phi_med = NA_real_,
             err = "")
}

append_row <- function(row) {
  write.table(row, CSV_OUT, sep = ",", row.names = FALSE,
              col.names = !file.exists(CSV_OUT), append = file.exists(CSV_OUT))
}

## ---- main loop ---------------------------------------------------------------
if (file.exists(CSV_OUT)) file.remove(CSV_OUT)

base_fits  <- list(A = NULL, B = NULL, C = NULL)   # compiled once, reused
data_hash  <- character(R)                          # cross-rep distinctness

t_start <- Sys.time()
for (r in seq_len(R)) {
  seed <- seeds[r]
  sim  <- sim_cohort(seed)
  data_hash[r] <- paste(round(sort(sim$t_dead), 10), collapse = "|")

  dat_naive <- make_naive(sim)
  dat_pp    <- make_pp(sim)

  cat(sprintf("[rep %02d seed %d] deaths observed: %d / %d;  pp rows: %d\n",
              r, seed, sum(sim$t_dead <= check_days), nrow(sim), nrow(dat_pp)))

  specs <- list(
    A = list(formula = form_naive, data = dat_naive,
             family = binomial(link = "identity"),
             prior = priors_mean, init = init_mean, has_phi = FALSE),
    B = list(formula = form_naive, data = dat_naive,
             family = brms::beta_binomial(link = "identity"),
             prior = priors_bb, init = init_bb, has_phi = TRUE),
    C = list(formula = form_pp, data = dat_pp,
             family = bernoulli(link = "identity"),
             prior = priors_mean, init = init_mean, has_phi = FALSE)
  )

  for (m in names(specs)) {
    sp  <- specs[[m]]
    row <- row_template(r, seed, m)
    res <- tryCatch({
      fit <- if (is.null(base_fits[[m]])) {
        brm(sp$formula, data = sp$data, family = sp$family, prior = sp$prior,
            chains = 1, iter = 800, warmup = 400, init = sp$init,
            seed = seed, refresh = 0, backend = "cmdstanr", silent = 2)
      } else {
        update(base_fits[[m]], newdata = sp$data, recompile = FALSE,
               chains = 1, iter = 800, warmup = 400, init = sp$init,
               seed = seed, refresh = 0, silent = 2)
      }
      if (is.null(base_fits[[m]])) base_fits[[m]] <- fit
      met <- extract_metrics(fit, sp$has_phi)
      if (!all(is.finite(c(met$z_med, met$z_lo, met$z_hi)))) {
        stop("non-finite z posterior summaries")
      }
      met
    }, error = function(e) e)

    if (inherits(res, "error")) {
      row$err <- conditionMessage(res)
      cat(sprintf("  method %s: FAILED (%s)\n", m, row$err))
    } else {
      row$fit_ok <- 1L
      for (nm in names(res)) row[[nm]] <- res[[nm]]
      cat(sprintf(
        "  method %s: z_med %.3f  CrI [%.3f, %.3f]  width %.3f  covers %d  div %g  Rhat %.3f%s\n",
        m, row$z_med, row$z_lo, row$z_hi, row$width, row$covers,
        row$n_div, row$max_rhat,
        if (sp$has_phi) sprintf("  phi_med %.1f", row$phi_med) else ""))
    }
    append_row(row)
  }
}
stopifnot("simulated cohorts distinct across replicates" =
            !anyDuplicated(data_hash))
t_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

## ---- aggregate & report ------------------------------------------------------
res <- read.csv(CSV_OUT, stringsAsFactors = FALSE)
stopifnot(nrow(res) == R * 3L)

method_label <- c(A = "A naive-binomial (cumulative cells)",
                  B = "B naive-betabinomial (package default family)",
                  C = "C correct person-period Bernoulli")

agg <- do.call(rbind, lapply(c("A", "B", "C"), function(m) {
  d  <- res[res$method == m, ]
  ok <- d[d$fit_ok == 1L, ]
  data.frame(
    method       = m,
    n_ok         = nrow(ok),
    n_failed     = sum(d$fit_ok == 0L),
    mean_z_med   = mean(ok$z_med),
    sd_z_med     = sd(ok$z_med),
    mean_width   = mean(ok$width),
    coverage     = sum(ok$covers),
    reps_div     = sum(ok$n_div > 0),
    total_div    = sum(ok$n_div),
    n_rhat_gt105 = sum(ok$max_rhat > 1.05),
    worst_rhat   = max(ok$max_rhat),
    min_ess      = min(ok$min_ess),
    b1_nonneg    = sum(ok$n_b1_nonneg)
  )
}))

## width ratios vs C: ratio of mean widths AND mean of per-rep ratios
wC <- res$width[res$method == "C"]
ratio_tab <- do.call(rbind, lapply(c("A", "B"), function(m) {
  wm   <- res$width[res$method == m]
  keep <- res$fit_ok[res$method == m] == 1L & res$fit_ok[res$method == "C"] == 1L
  data.frame(method = m,
             ratio_of_means = mean(wm[keep]) / mean(wC[keep]),
             mean_per_rep_ratio = mean(wm[keep] / wC[keep]),
             range_lo = min(wm[keep] / wC[keep]),
             range_hi = max(wm[keep] / wC[keep]))
}))

phi_ok <- res[res$method == "B" & res$fit_ok == 1L, "phi_med"]

old <- options(width = 110)
sink(TXT_OUT)
cat("==============================================================================\n")
cat("naive_vs_conditional_sim: cumulative-cell reformat vs correct person-period\n")
cat("likelihood for repeated-measures (same-cohort) survival checks\n")
cat("Run date:", format(Sys.time()), " |  total wall time:",
    sprintf("%.1f min", t_total), "\n")
cat("==============================================================================\n\n")

cat("--- Design (per replicate) ---------------------------------------------------\n")
cat("Temps 37/39/41/43 C (centred at 40), 40 individuals per temp, daily checks\n")
cat("days 1..14, right-censored at day 14. Truth: log10 t50 = 0.5 - 0.25*temp_c\n")
cat("(z = 4 C), log-logistic k = 8. R =", R, "replicates, seeds",
    min(seeds), "-", max(seeds), "(seed drives data AND sampler).\n")
cat("All fits: same 4PL mean structure (logk ~ 1, mid ~ temp_c), same priors\n")
cat("logk~N(1.5,1), mid_b0~N(0.5,1), mid_b1~N(0,0.5); B adds phi~gamma(2,0.1)\n")
cat("(package default). 1 chain, iter 800 (400 warmup), cmdstanr, no caching;\n")
cat("each model compiled once and reused via update(newdata=).\n\n")

cat("--- Methods ------------------------------------------------------------------\n")
cat("A: 56 rows temp x checkpoint, alive_j | trials(40), binomial(identity),\n")
cat("   mean S(t_j) -- cumulative counts treated as independent cells\n")
cat("B: same rows, beta_binomial(identity), phi estimated\n")
cat("C: person-period conditional Bernoulli, p_j = 1 - S(t_j)/S(t_{j-1})\n")
cat("   (log1p_exp spelling; algebraically = interval-censored log-logistic\n")
cat("   likelihood, tte_math_check.R claim (b))\n\n")

cat("--- Per-method summary (over fits that returned; failures counted, not dropped)\n")
for (i in seq_len(nrow(agg))) {
  a <- agg[i, ]
  cat(sprintf("%s\n", method_label[[a$method]]))
  cat(sprintf("  ok fits            : %d / %d   (failed: %d)\n", a$n_ok, R, a$n_failed))
  cat(sprintf("  mean z_med         : %.3f  (truth 4; sd across reps %.3f)\n",
              a$mean_z_med, a$sd_z_med))
  cat(sprintf("  mean 95%% CrI width : %.3f\n", a$mean_width))
  cat(sprintf("  coverage           : %d / %d\n", a$coverage, R))
  cat(sprintf("  divergences        : %d rep(s) with >0 (total %d)\n", a$reps_div, a$total_div))
  cat(sprintf("  max split-Rhat     : worst %.4f; reps > 1.05: %d; min bulk-ESS %.0f\n",
              a$worst_rhat, a$n_rhat_gt105, a$min_ess))
  cat(sprintf("  b1 >= 0 draws      : %d (z = -1/b1 sign-defined in all draws iff 0)\n",
              a$b1_nonneg))
  cat("\n")
}

cat("--- CrI width ratios vs method C ---------------------------------------------\n")
for (i in seq_len(nrow(ratio_tab))) {
  rt <- ratio_tab[i, ]
  cat(sprintf("%s / C: ratio of mean widths %.3f; mean per-rep ratio %.3f (range %.3f-%.3f)\n",
              rt$method, rt$ratio_of_means, rt$mean_per_rep_ratio,
              rt$range_lo, rt$range_hi))
}
cat("\n")

cat("--- Beta-binomial phi (method B) ---------------------------------------------\n")
cat(sprintf("posterior-median phi across reps: median %.1f, range %.1f-%.1f\n",
            median(phi_ok), min(phi_ok), max(phi_ok)))
cat("(large phi = beta-binomial collapsing toward plain binomial: within-cell\n")
cat(" overdispersion cannot represent BETWEEN-checkpoint serial dependence --\n")
cat(" marginally each cumulative count is exactly Binomial(40, S(t_j)).)\n\n")

cat("--- Per-replicate results ----------------------------------------------------\n")
print(res, row.names = FALSE)
cat("\n")

cat("--- What a 20-replicate study cannot rule out --------------------------------\n")
cat("* Coverage granularity: with R = 20 the MC standard error on a true 95%\n")
cat("  coverage is ~5 percentage points; small departures from nominal are not\n")
cat("  resolvable, only gross under-coverage.\n")
cat("* Single design point: one (z, k, b0), one temperature set, one n (40), one\n")
cat("  check schedule (daily, censor day 14). Overconfidence ratios can differ\n")
cat("  with sparser checks, other k, or heavier censoring.\n")
cat("* 1 chain x 400 post-warmup draws: within-chain split-Rhat only; CrI\n")
cat("  endpoints carry non-trivial MC error at 400 draws.\n")
cat("* Priors weakly informative and centred near truth; prior sensitivity and\n")
cat("  prior-data conflict not assessed.\n")
cat("* Generator matches the fitted mean model exactly (no frailty, no tank\n")
cat("  effects, no model misspecification); results isolate the dependence\n")
cat("  error alone.\n")
cat("* z-scale CrIs obtained by transforming b1 draws (z = -1/b1); with other\n")
cat("  designs where b1 posteriors approach 0 the transform itself misbehaves\n")
cat("  (b1 >= 0 draw counts reported above).\n")
sink()
options(old)

cat("\nDone. CSV:", CSV_OUT, "\nSummary:", TXT_OUT, "\n")
cat(sprintf("Total wall time: %.1f min\n", t_total))
