#!/usr/bin/env Rscript
# ==============================================================================
# tte_absolute_check.R
#
# Question: does bayesTLS's ABSOLUTE-threshold machinery (tls(target_surv=))
# work UNCHANGED on a person-period time-to-event fit, when the asymptotes are
# deliberately != 0/1 so relative and absolute thresholds genuinely differ?
#
# Design (constant-T cohorts, as in tte_prototype.R but with asymptotes):
#   4 temps (37, 39, 41, 43 C), n = 60 / temp, daily checks, censor day 14.
#   True 4PL survival-vs-log10(t):
#     S(t; T) = low + (up - low) / (1 + exp(k * (log10 t - mid(T))))
#     low = 0.05  (5% never die -> censored / cure fraction)
#     up  = 0.95  (5% die before the first check -> death in interval 1)
#     k   = 8,  mid(T) = 0.5 - 0.25 * temp_c  (temp_c = T - 40; true rel z = 4)
#   Because low/up/k are constant in T, the true ABSOLUTE z (any p) is also 4
#   and the absolute-vs-relative log10-LT offset is constant in T. NOTE:
#   low + up = 1 (symmetric about 0.5), so at p = 0.5 the offset is EXACTLY 0;
#   we therefore also read the curve at target_surv = 0.25, where the offset is
#   nonzero. All truths below are computed NUMERICALLY from the true S(t)
#   (uniroot), never from the closed-form correction (which is only used as an
#   independent cross-check).
#
# Fits (cmdstanr, 1 chain, iter 1000, NO file caching):
#   Fit A: FULL package-style nonlinear formula, all four nlpars
#          (lowraw, upraw, logk, mid), replicating make_4pl_formula()'s
#          disjoint-bounds inv_logit transform (compute_4pl_bounds via
#          getFromNamespace) and make_4pl_priors()'s midpoint priors.
#          Likelihood: person-period bernoulli(identity),
#          p_ij = 1 - S(lt_hi)/S_denom, S_denom = 1 for the FIRST interval
#          (conditions on S(0-) = 1, so the up-atom is observable in
#          interval 1) and S(lt_lo) otherwise.
#   Fit B: MINIMAL 2-nlpar recipe (logk, mid; low/up hard-coded 0/1).
#
# Then bayesTLS::tls() on each fit: relative, "absolute" (p = .5), p = .25.
# Plus derive_tdt_curve() smoke on both thresholds via a minimal bayes_tls
# wrapper (pattern: structure(list(fit=, data=, meta=), class="bayes_tls")).
#
# Scratchpad-only. Single fits: NO coverage/calibration claims -- recovery
# gates are deliberately loose existence checks; machinery gates (draw-level
# identities) are tight.
# ==============================================================================

SCRATCH <- "/private/tmp/claude-501/-Users-noble-Library-CloudStorage-Dropbox-1-Research-1-Manuscripts-In-Preparation-tls-model-equivalence/5004b512-34de-4b37-8d53-7208be821650/scratchpad"
SUMMARY_TXT <- file.path(SCRATCH, "tte_absolute_check_summary.txt")

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
  library(bayesTLS)
})
options(brms.backend = "cmdstanr", width = 100)

set.seed(20260926)
say <- function(...) cat(sprintf(...), "\n")

## Condition collector: capture value + messages + warnings + error text.
collect <- function(code) {
  msgs <- character(); warns <- character(); err <- NULL
  val <- withCallingHandlers(
    tryCatch(code, error = function(e) { err <<- conditionMessage(e); NULL }),
    message = function(m) { msgs  <<- c(msgs,  conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  list(value = val, messages = msgs, warnings = warns, error = err)
}
cond_report <- function(x, label) {
  say("  [%s] %s", label,
      if (is.null(x$error)) "SUCCEEDED" else paste0("ERROR: ", trimws(x$error)))
  for (m in x$messages) say("    message: %s", trimws(m))
  for (w in x$warnings) say("    warning: %s", trimws(w))
}

sink(SUMMARY_TXT, split = TRUE)
say("==============================================================================")
say("tte_absolute_check: absolute vs relative thresholds on a person-period TTE fit")
say("Run: %s | bayesTLS %s | brms %s", format(Sys.time()),
    as.character(utils::packageVersion("bayesTLS")),
    as.character(utils::packageVersion("brms")))
say("==============================================================================")

## =============================================================================
## 1. Truth, computed NUMERICALLY from the true S(t)
## =============================================================================
temps      <- c(37, 39, 41, 43)
t_mean     <- mean(temps)                      # 40 C centring constant
n_per_temp <- 60
check_days <- 14

true_low <- 0.05; true_up <- 0.95; true_k <- 8
true_b0  <- 0.5;  true_b1 <- -0.25            # mid(T) = b0 + b1 * temp_c
true_z   <- 4

S_true <- function(t, temp) {                  # t in days, temp in C
  mid <- true_b0 + true_b1 * (temp - t_mean)
  true_low + (true_up - true_low) / (1 + exp(true_k * (log10(t) - mid)))
}

## log10 time at which S crosses p, at temperature T (numeric root, no formula)
lt_at_p <- function(p, temp) {
  stopifnot(p > true_low, p < true_up)
  uniroot(function(lt) S_true(10^lt, temp) - p,
          interval = c(-8, 8), tol = 1e-12)$root
}

## Relative threshold value = (low+up)/2 = 0.5 here (symmetric asymptotes)
rel_thresh <- (true_low + true_up) / 2

Tg <- seq(36, 44, by = 0.5)
lt_rel   <- vapply(Tg, function(T) lt_at_p(rel_thresh, T), 0)
lt_a50   <- vapply(Tg, function(T) lt_at_p(0.50,       T), 0)
lt_a25   <- vapply(Tg, function(T) lt_at_p(0.25,       T), 0)
mid_true <- true_b0 + true_b1 * (Tg - t_mean)

off50 <- lt_a50 - mid_true                     # absolute(p=.5) - relative offset
off25 <- lt_a25 - mid_true                     # absolute(p=.25) - relative offset

## z from the numeric curves: -1/slope (least-squares slope over Tg; the curves
## are linear in T so this equals the local slope everywhere)
z_num <- function(lt) -1 / unname(coef(lm(lt ~ Tg))[2])
z_rel_true <- z_num(lt_rel); z_a50_true <- z_num(lt_a50); z_a25_true <- z_num(lt_a25)

## CTmax truths at t_ref = 1 day (log10 t = 0), by numeric root in T
ct_at <- function(lt_fun) uniroot(function(T) lt_fun(T), c(30, 50), tol = 1e-12)$root
ct_rel_true <- ct_at(function(T) lt_at_p(rel_thresh, T))
ct_a50_true <- ct_at(function(T) lt_at_p(0.50, T))
ct_a25_true <- ct_at(function(T) lt_at_p(0.25, T))

## Independent closed-form cross-check (NOT the source of truth):
off25_closed <- log((true_up - 0.25) / (0.25 - true_low)) / true_k

say("")
say("--- 1. Numeric truth from S(t) (t_ref = 1 day) -------------------------------")
say("relative threshold value (low+up)/2 = %.3f", rel_thresh)
say("offset  abs(p=.50) - rel : %.10f  (constant over T: sd %.2e)",
    mean(off50), sd(off50))
say("offset  abs(p=.25) - rel : %.10f  (constant over T: sd %.2e)",
    mean(off25), sd(off25))
say("closed-form ln((up-.25)/(.25-low))/k = %.10f  |diff vs numeric| = %.2e",
    off25_closed, abs(mean(off25) - off25_closed))
say("z: relative %.8f | abs p=.50 %.8f | abs p=.25 %.8f  (true z = 4)",
    z_rel_true, z_a50_true, z_a25_true)
say("CTmax(1 d): relative %.6f | abs p=.50 %.6f | abs p=.25 %.6f",
    ct_rel_true, ct_a50_true, ct_a25_true)

stopifnot(
  "truth: abs(p=.5) offset is exactly 0 (symmetric asymptotes)" =
    max(abs(off50)) < 1e-8,
  "truth: abs(p=.25) offset is constant in T" = sd(off25) < 1e-8,
  "truth: numeric abs(p=.25) offset matches closed form" =
    abs(mean(off25) - off25_closed) < 1e-8,
  "truth: absolute z = relative z = 4 (low/up/k constant in T)" =
    max(abs(c(z_rel_true, z_a50_true, z_a25_true) - true_z)) < 1e-6,
  "truth: relative CTmax(1d) = 42" = abs(ct_rel_true - 42) < 1e-6,
  "truth: abs(p=.25) CTmax = 42 + z*offset" =
    abs(ct_a25_true - (42 + true_z * mean(off25))) < 1e-6
)
say("TRUTH GATES PASSED")

## =============================================================================
## 2. Simulate cohort + person-period rows
## =============================================================================
## Inverse draw: U ~ unif(0,1) is the survival level at death.
##   U >= up  -> death at t = 0 (atom, mass 1-up; found dead at first check)
##   U <= low -> immortal (cure fraction, mass low; right-censored)
##   else     -> S(t) = U  =>  log10 t = mid + ln((up-U)/(U-low))/k
sim <- do.call(rbind, lapply(seq_along(temps), function(i) {
  Tt     <- temps[i]
  temp_c <- Tt - t_mean
  mid    <- true_b0 + true_b1 * temp_c
  U      <- runif(n_per_temp)
  ## branch explicitly (a vectorised ifelse would evaluate the log for the
  ## atom/immortal draws too and emit spurious NaN warnings)
  t_dead <- numeric(n_per_temp)
  t_dead[U >= true_up]  <- 0                                   # atom at t = 0
  t_dead[U <= true_low] <- Inf                                 # immortal (cure)
  mid_i <- U > true_low & U < true_up
  t_dead[mid_i] <- 10^(mid + log((true_up - U[mid_i]) /
                                   (U[mid_i] - true_low)) / true_k)
  data.frame(id = sprintf("T%d_%02d", i, seq_len(n_per_temp)),
             temp = Tt, temp_c = temp_c, t_dead = t_dead)
}))
sim$death_int <- ifelse(sim$t_dead <= check_days,
                        pmax(ceiling(sim$t_dead), 1L), NA)   # t=0 atom -> interval 1
sim$censored  <- is.na(sim$death_int)

## Person-period rows. FIRST interval conditions on S(0-) = 1: first = 1 and the
## denominator in the model is (first + (1-first)*S(lt_lo)), so lt_lo for j = 1
## is a finite placeholder (0) that is multiplied by zero -- NEVER a -Inf.
pp <- do.call(rbind, lapply(seq_len(nrow(sim)), function(r) {
  x <- sim[r, ]
  n_int <- if (x$censored) check_days else x$death_int
  j <- seq_len(n_int)
  data.frame(id = x$id, temp = x$temp, temp_c = x$temp_c,
             t_lo = j - 1, t_hi = j,
             lt_lo = ifelse(j == 1, 0, log10(j - 1)),   # placeholder when first=1
             lt_hi = log10(j),
             first = as.integer(j == 1),
             death = as.integer(!x$censored & j == n_int))
}))

say("")
say("--- 2. Simulated data --------------------------------------------------------")
say("individuals %d | person-period rows %d | deaths %d | censored %d | atom (t=0) deaths %d",
    nrow(sim), nrow(pp), sum(pp$death), sum(sim$censored), sum(sim$t_dead == 0))
print(with(sim, table(temp, censored)))
say("deaths in interval 1 by temp (contains the 1-up atom):")
print(with(sim[!sim$censored & sim$death_int == 1, ], table(temp)))

## =============================================================================
## 3. Fit A: FULL package-style 4-nlpar formula
## =============================================================================
## Replicate make_4pl_formula()'s disjoint-bounds transform exactly (same %.6f
## formatting), bounds (0, 1):
cb <- getFromNamespace("compute_4pl_bounds", "bayesTLS")
b  <- cb(0, 1)
low_expr <- sprintf("(%.6f + inv_logit(lowraw) * %.6f)", b$low_min, b$low_w)
up_expr  <- sprintf("(%.6f + inv_logit(upraw)  * %.6f)", b$up_min,  b$up_w)
S_expr   <- function(lt) sprintf(
  "(%s + (%s - %s) / (1 + exp(exp(logk) * (%s - mid))))",
  low_expr, up_expr, low_expr, lt)

rhs_A <- sprintf("1 - %s / (first + (1 - first) * %s)",
                 S_expr("lt_hi"), S_expr("lt_lo"))
## Clamped fallback (same quantity, guarded against p underflowing to exactly
## 0/1 during warmup; fmin/fmax are Stan builtins):
rhs_A_clamped <- sprintf("fmin(fmax(%s, 1e-12), 1 - 1e-12)", rhs_A)

form_A <- bf(as.formula(paste("death ~", rhs_A)),
             lowraw ~ 1, upraw ~ 1, logk ~ 1, mid ~ temp_c, nl = TRUE)
form_A_clamped <- bf(as.formula(paste("death ~", rhs_A_clamped)),
                     lowraw ~ 1, upraw ~ 1, logk ~ 1, mid ~ temp_c, nl = TRUE)

## Priors: make_4pl_priors() midpoint-branch centres/scales (bounds (0,1)):
##   lowraw ~ N(qlogis((0.02 - low_min)/low_w), 1); upraw likewise at 0.98;
##   logk ~ N(log 2, 1); mid Intercept ~ N(median(logd), 1.5); slopes N(0, 0.6).
## (No phi: bernoulli. No general class-b catch-alls: every coef gets an
##  explicit prior here, matching the coefficients that exist.)
lowraw_mean <- qlogis((0.02 - b$low_min) / b$low_w)
upraw_mean  <- qlogis((0.98 - b$up_min)  / b$up_w)
mid_start   <- median(pp$lt_hi)
priors_A <- c(
  set_prior(sprintf("normal(%.6f, 1)", lowraw_mean), class = "b",
            nlpar = "lowraw", coef = "Intercept"),
  set_prior(sprintf("normal(%.6f, 1)", upraw_mean),  class = "b",
            nlpar = "upraw",  coef = "Intercept"),
  set_prior(sprintf("normal(%.6f, 1)", log(2)),      class = "b",
            nlpar = "logk",   coef = "Intercept"),
  set_prior(sprintf("normal(%.6f, 1.5)", mid_start), class = "b",
            nlpar = "mid",    coef = "Intercept"),
  set_prior("normal(0, 0.6)", class = "b", nlpar = "mid", coef = "temp_c")
)

init_A <- function() list(b_lowraw = as.array(lowraw_mean),
                          b_upraw  = as.array(upraw_mean),
                          b_logk   = as.array(1.5),
                          b_mid    = c(0.5, 0))

fit_brm <- function(formula, priors, init, label) {
  message("Fitting: ", label)
  brm(formula, data = pp, family = bernoulli(link = "identity"),
      prior = priors, chains = 1, iter = 1000, warmup = 500,
      init = init, seed = 20260926, refresh = 0, backend = "cmdstanr")
}

say("")
say("--- 3. Fit A: 4-nlpar package-style person-period fit ------------------------")
say("main formula RHS: %s", rhs_A)
route_A <- "naive ratio spelling"
fitA <- tryCatch(fit_brm(form_A, priors_A, init_A, "Fit A (naive spelling)"),
                 error = function(e) e)
if (inherits(fitA, "error") || any(!is.finite(fixef(fitA)[, "Estimate"]))) {
  say("naive spelling failed (%s); refitting with fmin/fmax clamp",
      if (inherits(fitA, "error")) conditionMessage(fitA) else "non-finite estimates")
  route_A <- "fmin/fmax-clamped spelling of the same p"
  fitA <- fit_brm(form_A_clamped, priors_A, init_A, "Fit A (clamped spelling)")
}

dA <- as_draws_df(fitA)
low_d <- b$low_min + plogis(dA$b_lowraw_Intercept) * b$low_w
up_d  <- b$up_min  + plogis(dA$b_upraw_Intercept)  * b$up_w
k_d   <- exp(dA$b_logk_Intercept)
z_d   <- -1 / dA$b_mid_temp_c

npar <- nuts_params(fitA)
n_div <- sum(npar$Value[npar$Parameter == "divergent__"])
st <- summarise_draws(subset_draws(dA, variable = "^b_", regex = TRUE))
say("route used: %s | divergences %d | max Rhat %.4f | min bulk-ESS %.0f",
    route_A, n_div, max(st$rhat, na.rm = TRUE), min(st$ess_bulk, na.rm = TRUE))
print(st)
say("recovered (posterior medians): low %.4f (true .05) | up %.4f (true .95) | k %.3f (true 8) | z %.3f (true 4)",
    median(low_d), median(up_d), median(k_d), median(z_d))

stopifnot(
  "Fit A: asymptotes recovered to loose single-fit tolerance" =
    median(low_d) > 0.005 && median(low_d) < 0.2 &&
    median(up_d)  > 0.8   && median(up_d)  < 0.995,
  "Fit A: z within 25% of truth (single-fit existence gate)" =
    abs(median(z_d) - true_z) / true_z < 0.25
)

## =============================================================================
## 4. tls() on Fit A, both invocation styles, three thresholds
## =============================================================================
## Durations are DAYS and the fit is a bare brmsfit (meta = list()), so tls()'s
## t_ref/time_multiplier convention is: time_multiplier = 1 keeps t_ref in the
## MODEL's own units (days). We pass t_ref = 1, time_multiplier = 1 => CTmax at
## a 1-DAY reference exposure (log10 t_ref = 0), matching the truths above.
say("")
say("--- 4. tls() on Fit A (bare brmsfit) -----------------------------------------")
say("invocation: tls(fitA, temp_mean = 40, t_ref = 1, time_multiplier = 1,")
say("                temp_grid = seq(-4, 6, by = 0.5),   # extended: absolute CTmax")
say("                params = c('z','ctmax'), target_surv = <'relative'|'absolute'|0.25>)")
say("(default temp_grid spans only the observed temp_c range -3..3 = 37-43 C;")
say(" the absolute-p=.25 CTmax sits near 42.6 C, so draws crossing above 43 C")
say(" would fall off the default grid -> NA. tls()'s own docs say to extend it.)")

tgrid <- seq(-4, 6, by = 0.5)
tls_rel <- collect(tls(fitA, temp_mean = t_mean, t_ref = 1, time_multiplier = 1,
                       temp_grid = tgrid, params = c("z", "ctmax")))
tls_a50 <- collect(tls(fitA, temp_mean = t_mean, t_ref = 1, time_multiplier = 1,
                       temp_grid = tgrid, params = c("z", "ctmax"),
                       target_surv = "absolute"))
tls_a25 <- collect(tls(fitA, temp_mean = t_mean, t_ref = 1, time_multiplier = 1,
                       temp_grid = tgrid, params = c("z", "ctmax"),
                       target_surv = 0.25))
cond_report(tls_rel, "tls relative")
cond_report(tls_a50, "tls absolute p=0.50")
cond_report(tls_a25, "tls absolute p=0.25")
stopifnot("tls() ran on the person-period fit for all three thresholds" =
            !is.null(tls_rel$value) && !is.null(tls_a50$value) && !is.null(tls_a25$value))

get_med <- function(res, q) {
  s <- res$value$summary
  unlist(s[s$quantity == q, c("median", "lower", "upper")])
}
get_draws <- function(res, q) {
  d <- res$value$draws
  d$value[d$quantity == q][order(d$.draw[d$quantity == q])]
}

show3 <- function(label, v, truth) {
  say("%-22s median %8.4f  95%% CrI [%8.4f, %8.4f]   truth %8.4f", label,
      v["median"], v["lower"], v["upper"], truth)
}
say("")
show3("z   relative",     get_med(tls_rel, "z"),     true_z)
show3("z   absolute .50",  get_med(tls_a50, "z"),     true_z)
show3("z   absolute .25",  get_med(tls_a25, "z"),     true_z)
show3("CTmax relative",    get_med(tls_rel, "CTmax"), ct_rel_true)
show3("CTmax absolute .50", get_med(tls_a50, "CTmax"), ct_a50_true)
show3("CTmax absolute .25", get_med(tls_a25, "CTmax"), ct_a25_true)

zr <- get_draws(tls_rel, "z");     za50 <- get_draws(tls_a50, "z")
za25 <- get_draws(tls_a25, "z")
cr <- get_draws(tls_rel, "CTmax"); ca50 <- get_draws(tls_a50, "CTmax")
ca25 <- get_draws(tls_a25, "CTmax")

## Draw-level machinery identities (tight gates -- these test the EXTRACTION,
## not the single fit): with low/up/k constant in T the fitted absolute logLT
## curve is mid(T) + const per draw, so
##   (i)  z draws are identical across thresholds,
##   (ii) CTmax_abs - CTmax_rel = z * ln((up-p)/(p-low))/k  per draw.
off25_d <- log((up_d - 0.25) / (0.25 - low_d)) / k_d
off50_d <- log((up_d - 0.50) / (0.50 - low_d)) / k_d
stopifnot("all z / CTmax draws finite (no grid-inversion NAs)" =
            all(is.finite(c(zr, za50, za25, cr, ca50, ca25))))
say("")
say("draw-level identities:")
say("  max |z_abs25 - z_rel|                      = %.3e", max(abs(za25 - zr)))
say("  max |z_abs50 - z_rel|                      = %.3e", max(abs(za50 - zr)))
say("  max |(CT_abs25 - CT_rel) - z*offset25(dr)| = %.3e",
    max(abs((ca25 - cr) - z_d * off25_d)))
say("  max |(CT_abs50 - CT_rel) - z*offset50(dr)| = %.3e",
    max(abs((ca50 - cr) - z_d * off50_d)))
say("  median fitted offset25 %.5f (true %.5f) | median fitted offset50 %.5f (true 0)",
    median(off25_d), mean(off25), median(off50_d))
say("  median z draws: rel %.6f = -1/b_mid_temp_c %.6f", median(zr), median(z_d))

stopifnot(
  "machinery: absolute z draws == relative z draws (offset constant in T)" =
    max(abs(za25 - zr)) < 1e-6 && max(abs(za50 - zr)) < 1e-6,
  "machinery: pooled relative z == -1/b_mid_temp_c per draw" =
    max(abs(zr - z_d)) < 1e-6,
  "machinery: CTmax_abs - CTmax_rel == z * ln((up-p)/(p-low))/k per draw" =
    max(abs((ca25 - cr) - z_d * off25_d)) < 1e-6 &&
    max(abs((ca50 - cr) - z_d * off50_d)) < 1e-6,
  "recovery: absolute-.25 CTmax median within 1 C of truth (loose, single fit)" =
    abs(unname(get_med(tls_a25, "CTmax")["median"]) - ct_a25_true) < 1,
  "recovery: relative CTmax median within 1 C of truth (loose, single fit)" =
    abs(unname(get_med(tls_rel, "CTmax")["median"]) - ct_rel_true) < 1
)
say("FIT-A MACHINERY GATES PASSED")

## ----- workflow-wrapper invocation (minimal bayes_tls wrapper) ---------------
## Wrapper meta gives tls() everything it otherwise needs as arguments:
## temp_mean, duration_unit = days (-> time_multiplier 1440 to minutes),
## t_ref recorded in MINUTES (1440 = 1 day).
wrapA <- structure(list(
  fit = fitA, data = transform(pp, duration = t_hi),
  meta = list(temp_mean = t_mean, duration_unit = "days", t_ref = 1440,
              response_type = "proportion", parameterization = "midpoint",
              group_vars = character(0), grouped = FALSE,
              random_effects = NULL, lower = 0, upper = 1,
              threshold = "relative", log10_tref = 0)),
  class = "bayes_tls")

tls_wrap <- collect(tls(wrapA, params = c("z", "ctmax"), target_surv = 0.25,
                        temp_grid = tgrid))
say("")
say("wrapper invocation: tls(wrapA, params = c('z','ctmax'), target_surv = 0.25,")
say("                        temp_grid = seq(-4, 6, by = 0.5))")
cond_report(tls_wrap, "tls on bayes_tls wrapper, p=0.25")
if (!is.null(tls_wrap$value)) {
  zw <- get_draws(tls_wrap, "z"); cw <- get_draws(tls_wrap, "CTmax")
  say("  wrapper vs bare-brmsfit draws: max|z diff| %.3e  max|CTmax diff| %.3e",
      max(abs(zw - za25)), max(abs(cw - ca25)))
  stopifnot("wrapper and bare-brmsfit tls() agree draw-for-draw" =
              max(abs(zw - za25)) < 1e-9 && max(abs(cw - ca25)) < 1e-9)
}

## =============================================================================
## 5. derive_tdt_curve() smoke on Fit A (wrapper; both thresholds)
## =============================================================================
say("")
say("--- 5. derive_tdt_curve() smoke (wrapper, temp_grid = c(38, 40, 42)) ---------")
dc_rel <- collect(derive_tdt_curve(wrapA, temp_grid = c(38, 40, 42)))
dc_a25 <- collect(derive_tdt_curve(wrapA, temp_grid = c(38, 40, 42),
                                   target_surv = 0.25))
cond_report(dc_rel, "derive_tdt_curve relative")
cond_report(dc_a25, "derive_tdt_curve absolute p=0.25")
if (!is.null(dc_rel$value)) print(dc_rel$value$summary)
if (!is.null(dc_a25$value)) {
  print(dc_a25$value$summary)
  say("  NOTE: even where it RUNS, the absolute branch predicts the model's mu,")
  say("  which for a person-period fit is the CONDITIONAL interval death")
  say("  probability 1 - S(hi)/S(lo), NOT survival S(t): the inverted 'LT' is")
  say("  then the wrong quantity. Only the nlpar-based readers (tls) are")
  say("  semantically safe on a person-period fit without changes.")
}

## =============================================================================
## 6. Fit B: minimal 2-nlpar recipe (low/up hard-coded 0/1)
## =============================================================================
say("")
say("--- 6. Fit B: minimal 2-nlpar recipe (logk, mid; low = 0, up = 1) ------------")
SB <- function(lt) sprintf("(1 / (1 + exp(exp(logk) * (%s - mid))))", lt)
rhs_B <- sprintf("1 - %s / (first + (1 - first) * %s)", SB("lt_hi"), SB("lt_lo"))
form_B <- bf(as.formula(paste("death ~", rhs_B)),
             logk ~ 1, mid ~ temp_c, nl = TRUE)
priors_B <- c(
  set_prior(sprintf("normal(%.6f, 1)", log(2)),      class = "b",
            nlpar = "logk", coef = "Intercept"),
  set_prior(sprintf("normal(%.6f, 1.5)", mid_start), class = "b",
            nlpar = "mid",  coef = "Intercept"),
  set_prior("normal(0, 0.6)", class = "b", nlpar = "mid", coef = "temp_c")
)
init_B <- function() list(b_logk = as.array(1.5), b_mid = c(0.5, 0))
say("main formula RHS: %s", rhs_B)
fitB <- fit_brm(form_B, priors_B, init_B, "Fit B (2-nlpar, asymptotes fixed 0/1)")
dB <- as_draws_df(fitB)
say("Fit B (misspecified for these data: ignores the 5%% atom + 5%% cure):")
say("  z = -1/b_mid_temp_c median %.3f (true 4) | k median %.3f (true 8)",
    median(-1 / dB$b_mid_temp_c), median(exp(dB$b_logk_Intercept)))

tlsB_rel <- collect(tls(fitB, temp_mean = t_mean, t_ref = 1, time_multiplier = 1,
                        temp_grid = tgrid, params = c("z", "ctmax")))
tlsB_a25 <- collect(tls(fitB, temp_mean = t_mean, t_ref = 1, time_multiplier = 1,
                        temp_grid = tgrid, params = c("z", "ctmax"),
                        target_surv = 0.25))
say("")
say("tls() on Fit B (default nlpars = c('lowraw','upraw','logk','mid')):")
cond_report(tlsB_rel, "tls relative on 2-nlpar fit")
cond_report(tlsB_a25, "tls absolute p=0.25 on 2-nlpar fit")
stopifnot("Fit B: tls() default errors (no lowraw/upraw nlpars in the fit)" =
            !is.null(tlsB_rel$error) && !is.null(tlsB_a25$error))

## Workaround probe: aliasing the missing asymptote nlpars onto an existing one
## makes RELATIVE mode run (logLT = mid; low/up/k unused downstream) -- but the
## same alias in ABSOLUTE mode would use garbage low/up, so it is a footgun,
## recorded here for the interface discussion, not a recommendation.
tlsB_alias_rel <- collect(tls(fitB, temp_mean = t_mean, t_ref = 1,
                              time_multiplier = 1, temp_grid = tgrid,
                              params = c("z", "ctmax"),
                              nlpars = c("logk", "logk", "logk", "mid")))
tlsB_alias_a25 <- collect(tls(fitB, temp_mean = t_mean, t_ref = 1,
                              time_multiplier = 1, temp_grid = tgrid,
                              params = c("z", "ctmax"), target_surv = 0.25,
                              nlpars = c("logk", "logk", "logk", "mid")))
say("")
say("alias probe: tls(fitB, nlpars = c('logk','logk','logk','mid')):")
cond_report(tlsB_alias_rel, "alias, relative")
if (!is.null(tlsB_alias_rel$value)) {
  show3("  z   (alias rel)",    get_med(tlsB_alias_rel, "z"),     true_z)
  show3("  CTmax (alias rel)",  get_med(tlsB_alias_rel, "CTmax"), ct_rel_true)
}
cond_report(tlsB_alias_a25, "alias, absolute p=0.25 (DO NOT TRUST values)")
if (!is.null(tlsB_alias_a25$value)) {
  say("  alias-absolute uses inv_logit(logk-draws) as FAKE asymptotes -> wrong")
  say("  threshold by construction. Whether it fails loudly (NA + warning, when")
  say("  p falls outside the fake asymptote range) or silently returns wrong")
  say("  numbers (when p happens to fall inside it) depends on the logk draws:")
  show3("  CTmax (alias abs .25)", get_med(tlsB_alias_a25, "CTmax"), ct_a25_true)
}

## =============================================================================
## 7. Verdict block
## =============================================================================
say("")
say("==============================================================================")
say("VERDICT")
say("==============================================================================")
say("* The RELATIVE/ABSOLUTE distinction lives entirely in post-fit extraction:")
say("  the fitted model parameter mid is the RELATIVE (midpoint) log10-LT; the")
say("  absolute LTp adds ln((up - p)/(p - low))/k per draw inside tls(), which")
say("  works UNCHANGED on the person-period TTE fit PROVIDED the fit carries all")
say("  four nlpars (Fit A). Same correction, same code path as the classical fit.")
say("* With low/up/k constant in T the absolute z equals the relative z exactly")
say("  (draw-for-draw); CTmax shifts by z * offset. Verified against numeric truth.")
say("* The minimal 2-nlpar recipe (Fit B) supports NEITHER tls() mode as-is:")
say("  posterior_linpred(nlpar='lowraw') fails. Relative-only extraction is")
say("  possible via the nlpars alias hack; absolute is NOT safely available")
say("  without asymptote nlpars in the fit.")
say("* Single simulated dataset, single fits, 1 chain: existence proof only;")
say("  no coverage/calibration claims.")
sink()
cat("\nSummary written to:", SUMMARY_TXT, "\n")
