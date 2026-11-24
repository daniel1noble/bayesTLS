# Shared definitions for the golden builder fixtures (test-golden-stancode.R).
# Kept in a helper so a baseline can be regenerated from another checkout of the
# package (load that package, source this file and helper-simulate.R, call
# golden_snapshot() on each config) without running the test itself.

golden_dir <- function() testthat::test_path("golden")

g_std <- function() {
  standardize_data(simulate_tdt(seed = 1), temp = "T", duration = "t_hours",
                   n_total = "n_trials", n_surv = "y_alive",
                   duration_unit = "hours")
}
g_std2 <- function() {
  standardize_data(simulate_tdt_2group(seed = 2), temp = "T",
                   duration = "t_hours", n_total = "n_trials",
                   n_surv = "y_alive", duration_unit = "hours")
}
g_beta <- function() {
  standardize_data(simulate_tdt_beta(seed = 3), temp = "T",
                   duration = "t_hours", proportion = "y_prop",
                   duration_unit = "hours")
}
g_std_batch <- function() {                # midpoint random intercept
  d <- simulate_tdt(seed = 1)
  d$batch <- factor(rep(1:5, length.out = nrow(d)))
  standardize_data(d, temp = "T", duration = "t_hours", n_total = "n_trials",
                   n_surv = "y_alive", duration_unit = "hours")
}
g_std_min <- function() {                  # minutes, as in every case study
  d <- simulate_tdt(seed = 1)
  d$t_min <- d$t_hours * 60
  standardize_data(d, temp = "T", duration = "t_min", n_total = "n_trials",
                   n_surv = "y_alive", duration_unit = "minutes")
}
g_beta_plant <- function() {               # snow-gum-like: Beta + direct + RE
  d <- simulate_tdt_beta_2group(seed = 4)
  d$plant <- factor(rep(1:6, length.out = nrow(d)))
  standardize_data(d, temp = "T", duration = "t_hours", proportion = "y_prop",
                   duration_unit = "hours")
}

# One entry per frozen configuration: the model it uses and how to build it.
# Log-logistic entries were generated from the package BEFORE the
# alternative-models work (commit 5a671dc), so they pin that behaviour; the
# other models' entries pin the first version that supported them.
golden_configs <- function() {
  ll <- function(f) list(model = "log-logistic", build = f)
  list(
    mid_bb_default = ll(function() fit_4pl(g_std(), fit = FALSE)),
    mid_bb_by      = ll(function() fit_4pl(g_std2(), by = "grp", fit = FALSE)),
    mid_bb_shape   = ll(function() fit_4pl(g_std(), low = ~ 1, up = ~ temp_c,
                                           k = ~ temp_c, mid = ~ temp_c,
                                           fit = FALSE)),
    mid_binomial   = ll(function() {
      fit_4pl(g_std(), family = stats::binomial(link = "identity"), fit = FALSE)
    }),
    mid_beta       = ll(function() fit_4pl(g_beta(), fit = FALSE)),
    direct_rel_1   = ll(function() fit_4pl(g_std(), ctmax = ~ 1, z = ~ 1,
                                           fit = FALSE)),
    direct_rel_grp = ll(function() fit_4pl(g_std2(), ctmax = ~ 0 + grp,
                                           z = ~ 0 + grp, fit = FALSE)),
    direct_abs_1   = ll(function() fit_4pl(g_std(), ctmax = ~ 1, z = ~ 1,
                                           threshold = "absolute", fit = FALSE)),
    mid_fixed_asy  = ll(function() fit_4pl(g_std(), fix_asymptotes = TRUE,
                                           fit = FALSE)),
    mid_re         = ll(function() fit_4pl(g_std_batch(),
                                           random_effects = "batch",
                                           fit = FALSE)),
    direct_re_beta = ll(function() fit_4pl(g_beta_plant(),
                                           ctmax = ~ 0 + grp + (1 | plant),
                                           z = ~ 0 + grp, fit = FALSE)),
    direct_min_tref240 = ll(function() fit_4pl(g_std_min(), ctmax = ~ 1,
                                               z = ~ 1, t_ref = 240,
                                               fit = FALSE)),
    mid_bounds_sub = ll(function() fit_4pl(g_std(), bounds = c(0, 0.9),
                                           fit = FALSE)),
    direct_abs_p045 = ll(function() fit_4pl(g_std(), ctmax = ~ 1, z = ~ 1,
                                            bounds = c(0, 0.9),
                                            threshold = "absolute", p = 0.45,
                                            fit = FALSE)),
    ln_direct_abs  = list(model = "log-normal", build = function() {
      fit_lognormal(g_std(), ctmax = ~ 1, z = ~ 1, threshold = "absolute",
                    fit = FALSE)
    }),
    w1_mid_default = list(model = "weibull1", build = function() {
      fit_weibull(g_std(), fit = FALSE)
    }),
    w2_direct_grp  = list(model = "weibull2", build = function() {
      fit_weibull(g_std2(), type = 2, ctmax = ~ 0 + grp, z = ~ 0 + grp,
                  fit = FALSE)
    })
  )
}

# brms-independent views: the formulas as base R deparses them (they come from
# our own strings via stats::as.formula()) and the prior table's defining
# columns, sorted. These are compared on any brms version.
golden_formula_parts <- function(bf) {
  dep <- function(f) paste(deparse(f, width.cutoff = 500L), collapse = " ")
  pf  <- bf$pforms[order(names(bf$pforms))]
  c(main = dep(bf$formula), vapply(pf, dep, character(1)),
    family = paste(bf$family$family, bf$family$link))
}
golden_prior_core <- function(pr) {
  d <- as.data.frame(pr)[, c("prior", "class", "coef", "group", "nlpar")]
  d <- d[do.call(order, unname(as.list(d))), , drop = FALSE]
  rownames(d) <- NULL
  d
}

golden_snapshot <- function(wf) {
  list(formula_parts = golden_formula_parts(wf$formula),
       prior_core    = golden_prior_core(wf$prior),
       # brms-rendered views: compared only on the brms version recorded here.
       # suppressWarnings: brms reports unused global `b` priors for
       # intercept-only direct configs (pre-existing, harmless).
       formula       = utils::capture.output(print(wf$formula)),
       prior         = as.data.frame(wf$prior),
       stancode      = suppressWarnings(
         brms::make_stancode(wf$formula, data = wf$data, prior = wf$prior,
                             data2 = wf$data2)),
       brms          = as.character(utils::packageVersion("brms")))
}
