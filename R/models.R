# Dose-response model registry.
#
# Every model is an S-shaped curve written in its MIDPOINT form, survival
# equal to low plus (up - low) times G of k (log10 t - mid),
# where `mid` is the log10 time at which p = (low + up)/2 for EVERY model, so
# the TDT line mid(T) = log10 t_ref + (CTmax - T)/z, the relative threshold,
# z and CTmax mean the same thing whichever curve is fitted. The entries hold
# the pieces the package needs to swap the curve without touching anything
# else:
#   rhs(low, up)          brms non-linear formula RHS (Stan syntax) given the
#                         expressions used for the two asymptotes
#   abs_corr_stan(p)      Stan expression for c/k, the correction that turns
#                         mid into the absolute log10 t at survival p (baked
#                         into direct fits with threshold = "absolute")
#   G(u)                  the curve on the R side, as a share of the range
#   curve(low, up, u)     the full curve, p at u = k (log10 t - mid); the
#                         log-logistic form is written exactly as the
#                         pre-models code wrote it, so its numbers are
#                         bit-identical
#   abs_corr(low,up,k,p)  the same correction on the R side
#   slope_factor          the slope of p against log10 t at mid equals
#                         minus k times (up - low) times slope_factor
#   theta_offset          standard (inflection-point) location of the curve
#                         is mid + theta_offset / k (reporting only)
#   hazard(t, p, k, mid)  hazard for low = 0, up = 1
#   signature             whitespace-free fragment of rhs() that no other
#                         model's formula contains; identifies the curve of a
#                         bare brms fit (see tdt_model_from_formula())
# The log-logistic entry reproduces the pre-1.2 strings character for
# character, so its Stan code, priors and cached fits are unchanged.

tdt_models <- function() {
  list(
    "log-logistic" = list(
      name          = "log-logistic",
      signature     = "/(1+exp(exp(logk)*(logd-mid)))",
      label         = "log-logistic (4PL)",
      aliases       = c("loglogistic", "4pl", "logistic", "logit"),
      rhs           = function(low, up) {
        sprintf("%s + (%s - %s) / (1 + exp(exp(logk) * (logd - mid)))",
                low, up, low)
      },
      abs_corr_stan = function(p) {
        sprintf("(1/exp(logk)) * log((up - %g)/(%g - low))", p, p)
      },
      G             = function(u) 1 / (1 + exp(u)),
      curve         = function(low, up, u) low + (up - low) / (1 + exp(u)),
      abs_corr      = function(low, up, k, p) log((up - p) / (p - low)) / k,
      slope_factor  = 0.25,
      theta_offset  = 0,
      hazard        = function(t, p, k, mid) k / (t * log(10)) * (1 - p),
      death_times   = "log-logistic"
    ),
    "log-normal" = list(
      name          = "log-normal",
      signature     = "*Phi(-exp(logk)*(logd-mid))",
      label         = "log-normal (probit)",
      aliases       = c("lognormal", "probit", "normal"),
      rhs           = function(low, up) {
        sprintf("%s + (%s - %s) * Phi(-exp(logk) * (logd - mid))",
                low, up, low)
      },
      abs_corr_stan = function(p) {
        sprintf("(1/exp(logk)) * inv_Phi((up - %g)/(up - low))", p)
      },
      G             = function(u) stats::pnorm(-u),
      curve         = function(low, up, u) low + (up - low) * stats::pnorm(-u),
      abs_corr      = function(low, up, k, p) {
        stats::qnorm((up - p) / (up - low)) / k
      },
      slope_factor  = stats::dnorm(0),
      theta_offset  = 0,
      hazard        = function(t, p, k, mid) {
        k / (t * log(10)) * stats::dnorm(k * (log10(t) - mid)) / p
      },
      death_times   = "log-normal"
    ),
    "weibull1" = list(
      name          = "weibull1",
      signature     = "*exp(-log(2)*exp(exp(logk)*(logd-mid)))",
      label         = "Weibull type 1 (Gompertz-shaped)",
      aliases       = c("weibull-1", "weibull_1", "w1", "weibull", "gompertz"),
      rhs           = function(low, up) {
        sprintf(paste0("%s + (%s - %s) * ",
                       "exp(-log(2) * exp(exp(logk) * (logd - mid)))"),
                low, up, low)
      },
      abs_corr_stan = function(p) {
        sprintf("(1/exp(logk)) * log(log2((up - low)/(%g - low)))", p)
      },
      G             = function(u) 2^(-exp(u)),
      curve         = function(low, up, u) low + (up - low) * 2^(-exp(u)),
      abs_corr      = function(low, up, k, p) {
        log(log2((up - low) / (p - low))) / k
      },
      slope_factor  = log(2) / 2,
      theta_offset  = -log(log(2)),
      hazard        = function(t, p, k, mid) k / (t * log(10)) * (-log(p)),
      death_times   = "Weibull"
    ),
    "weibull2" = list(
      name          = "weibull2",
      signature     = "*exp(-log(2)*exp(-exp(logk)*(logd-mid)))",
      label         = "Weibull type 2",
      aliases       = c("weibull-2", "weibull_2", "w2", "frechet"),
      rhs           = function(low, up) {
        sprintf(paste0("%s - (%s - %s) * ",
                       "exp(-log(2) * exp(-exp(logk) * (logd - mid)))"),
                up, up, low)
      },
      abs_corr_stan = function(p) {
        sprintf("(-(1/exp(logk)) * log(log2((up - low)/(up - %g))))", p)
      },
      G             = function(u) 1 - 2^(-exp(-u)),
      curve         = function(low, up, u) up - (up - low) * 2^(-exp(-u)),
      abs_corr      = function(low, up, k, p) {
        -log(log2((up - low) / (up - p))) / k
      },
      slope_factor  = log(2) / 2,
      theta_offset  = log(log(2)),
      hazard        = function(t, p, k, mid) {
        k / (t * log(10)) * (-(1 - p) * log(1 - p)) / p
      },
      death_times   = "Fr\u00e9chet (inverse Weibull)"
    )
  )
}

# Canonical model name from a user-facing name or alias (case, spaces, "-"
# and "_" ignored), or an error listing the choices.
tdt_match_model <- function(model) {
  if (is.null(model)) return("log-logistic")
  reg <- tdt_models()
  if (!is.character(model) || length(model) != 1L || is.na(model))
    stop("`model` must be one of: ", paste(names(reg), collapse = ", "), ".",
         call. = FALSE)
  norm <- function(x) tolower(gsub("[ _-]", "", x))
  for (nm in names(reg))
    if (norm(model) %in% norm(c(nm, reg[[nm]]$aliases))) return(nm)
  stop(sprintf("Unknown model \"%s\". Choose one of: %s.", model,
               paste(names(reg), collapse = ", ")), call. = FALSE)
}

tdt_model_spec <- function(model) tdt_models()[[tdt_match_model(model)]]

# Which model a brms formula encodes, from the curve fragment each registry
# entry carries; NA when none matches (a hand-written formula).
tdt_model_from_formula <- function(formula) {
  if (inherits(formula, "brmsformula")) formula <- formula$formula
  if (!inherits(formula, "formula")) return(NA_character_)
  txt <- gsub("[[:space:]]", "",
              paste(deparse(formula, width.cutoff = 500L), collapse = ""))
  reg <- tdt_models()
  hit <- names(reg)[vapply(reg, function(m) {
    grepl(m$signature, txt, fixed = TRUE)
  }, logical(1))]
  if (length(hit) == 1L) hit else NA_character_
}

# Which model an object was fitted with (canonical name). A workflow records it
# in `meta$model`; a workflow saved before models existed, or a bare brmsfit, is
# recognised from its formula. Anything unrecognised -- including a hand-written
# brms 4PL -- is treated as the log-logistic, the package's original model.
tdt_model <- function(x) {
  m <- NULL
  if (inherits(x, "bayes_tls")) {
    m <- x$meta$model
    if (is.null(m)) {
      f <- if (inherits(x$fit, "brmsfit")) x$fit$formula else x$formula
      if (!is.null(f)) m <- tdt_model_from_formula(f)
    }
  } else if (inherits(x, "brmsfit")) {
    m <- tdt_model_from_formula(x$formula)
  }
  if (is.null(m) || is.na(m)) "log-logistic" else tdt_match_model(m)
}

# A model spec from a name/alias or an already-resolved spec. Hot loops (the
# heat-injury integrators) resolve the spec once and pass the list along.
tdt_as_model_spec <- function(model) {
  if (is.list(model)) model else tdt_model_spec(model)
}

# Absolute offset, in log10 time, of survival target p from the midpoint:
# the absolute log10 LT at p is mid + offset. Element-wise over low, up, k
# (scalars, vectors or matrices of one shape; p a scalar). The target is
# reachable only strictly between the asymptotes, low < p < up; elsewhere the
# result is NA, and the model's correction is never evaluated there, so no
# NaN warnings arise. For the log-logistic this is log((up - p)/(p - low))/k,
# computed exactly as the pre-models engine computed it.
tdt_abs_offset <- function(model, low, up, k, p) {
  mod <- tdt_as_model_spec(model)
  n   <- max(length(low), length(up), length(k))
  lo  <- rep_len(low, n)
  hi  <- rep_len(up, n)
  kk  <- rep_len(k, n)
  ok  <- !is.na(lo) & !is.na(hi) & lo < p & p < hi
  out <- rep(NA_real_, n)
  out[ok] <- mod$abs_corr(lo[ok], hi[ok], kk[ok], p)
  dim(out) <- dim(low) %||% dim(up) %||% dim(k)
  out
}

# Prior mean for logk such that the implied prior on the slope of p at the
# midpoint, -k (up - low) * slope_factor, is the same for every model. The
# log-logistic keeps its original log(2) exactly (log(1) == 0).
tdt_logk_prior_mean <- function(model) {
  log(2) + log(0.25 / tdt_model_spec(model)$slope_factor)
}

# Inverse standard-normal CDF under its Stan name, used by the log-normal
# absolute-threshold correction. In Stan it is the built-in inv_Phi(). On the R
# side brms (prepare_predictions.bframenl) evaluates each non-linear formula in
# env_stan_functions(parent = environment(<formula>)), which supplies Phi and
# inv_logit but not inv_Phi or log2. The formula environment is the frame of
# make_tdt_formula() (stats::as.formula() is called there) and that frame's
# parent is this namespace, so inv_Phi resolves here and log2 in base. Keep the
# formulas built with stats::as.formula() inside make_tdt_formula(), and do not
# reset their environment.
inv_Phi <- function(x) stats::qnorm(x) # nolint: object_name_linter.
