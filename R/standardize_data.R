#' Standardise a raw survival / proportion dataset for the TDT function library
#'
#' Rewrites user column names into a single project-standard schema and attaches
#' metadata used by every downstream fitting and prediction helper. This is the
#' single entry point for raw data — everything else in the library assumes the
#' output of this function.
#'
#' Two response types are supported:
#'
#' - **Count data** (binomial / beta-binomial): supply `n_total` plus **exactly
#'   one** of `n_surv`, `n_dead`, `survival`, or `mortality`. The other counts
#'   are derived and the standardised columns include `n_total`, `n_surv`,
#'   `n_dead`, `survival`. `response_type` is recorded as `"count"`.
#' - **Continuous proportion** (Beta), e.g. a chlorophyll-fluorescence
#'   \eqn{F_v/F_m} ratio with no denominator: supply `proportion` and omit the
#'   count arguments. The value is stored in `survival` (clamped into the open
#'   interval `(proportion_eps, 1 - proportion_eps)` so the Beta likelihood is
#'   finite, except on a side with a `detection_limit`); no
#'   `n_total`/`n_surv` columns are created. `response_type` is
#'   recorded as `"proportion"`. If the instrument cannot resolve values
#'   beyond a limit (e.g. it records 0 when function is below what it can
#'   read), give that limit as `detection_limit` so those readings are treated
#'   as censored rather than clamped.
#'
#' If the dataset spans multiple categories (life stages, species, populations,
#' etc.), filter to one category before calling this function and fit a separate
#' model per subset — the fitter does not estimate category-level effects.
#'
#' @param data           Raw data frame or tibble.
#' @param temp           Column name of the assay temperature (°C).
#' @param duration       Column name of the exposure duration. The unit is
#'                       whatever is in the source data; record it via
#'                       `duration_unit`.
#' @param n_total        Column name for total individuals per replicate.
#'                       Required for count responses; leave `NULL` (default)
#'                       for a continuous `proportion` response.
#' @param n_surv         Column name for survivor counts.
#' @param n_dead         Column name for death counts. Converted to `n_surv`
#'                       via `n_surv = n_total - n_dead`.
#' @param survival       Column name for survival proportions in `[0, 1]`.
#'                       Converted to integer counts via `n_total`.
#' @param mortality      Column name for mortality proportions in `[0, 1]`.
#'                       Converted to `n_surv = round((1 - mortality) * n_total)`.
#' @param proportion     Column name for a continuous proportion response in
#'                       `[0, 1]` with no denominator (modelled with a Beta
#'                       likelihood). Mutually exclusive with the count
#'                       arguments above.
#' @param proportion_eps Half-open clamp applied to `proportion` so values sit
#'                       strictly inside `(0, 1)` (the Beta density is undefined
#'                       at exactly 0 or 1): one number in `(0, 0.5)`. Default
#'                       `0.001`. With a `detection_limit`, the clamp applies
#'                       only on a side that has no limit; censored readings
#'                       take their limit, and uncensored readings on a side
#'                       with a limit are kept as recorded.
#' @param detection_limit For a continuous `proportion` response only: the
#'   limit(s) beyond which the measurement cannot resolve values, so that a
#'   reading there says only that the true value lies beyond the limit (a
#'   censored reading). `NULL` (default): no censoring, and readings of 0 or 1
#'   are clamped by `proportion_eps`. Otherwise one of:
#'   - a number in `(0, 1)`: a lower limit. Readings below it (typically
#'     recorded as 0) are left-censored: the true value is below the limit.
#'   - a length-2 numeric `c(lower, upper)`, named or in that order, with `NA`
#'     for no limit on one side (or one named number, `c(upper = 0.98)`).
#'     Readings at or above `upper` (e.g. values capped at a ceiling) are
#'     right-censored: the true value is at or above the limit.
#'   - the name of a numeric column holding a per-row lower limit, or a
#'     length-2 character `c(lower = "col_lo", upper = "col_hi")` (a list
#'     mixes numbers and column names). `NA` in a limit column means no limit
#'     for that row.
#'
#'   Censored rows store their limit in `survival`, are flagged `"left"` or
#'   `"right"` in a new `censored` column (`"none"` otherwise), and enter the
#'   Beta likelihood as the probability of lying beyond the limit (brms
#'   `cens()`, added by [fit_4pl()] and the other fitters); the reading as
#'   recorded is kept in `survival_recorded` for plots and observed summaries.
#'   The fitted curve keeps its meaning (the mean of the uncensored response),
#'   so [tls()] and the prediction helpers are unchanged. There is no default
#'   limit: take it from the instrument or the data (e.g. the smallest positive
#'   reading). A reading equal to the lower limit counts as observed; if your
#'   instrument records below-limit readings AT the limit rather than as 0,
#'   recode those to 0 first. Positive readings below a lower limit, and
#'   readings strictly between an upper limit and 1, are coarsened to censored
#'   with a warning, since they suggest the limit is set too high (or too low).
#'
#'   The pointwise log-likelihood of a censored fit (`log_lik()`, `loo()`,
#'   `waic()`) scores each censored reading as a probability, not a density, so
#'   compare such fits only with fits to the same data and the same
#'   `detection_limit`, never with a fit whose readings were clamped.
#' @param random_effects Optional character vector of grouping variables for
#'                       random effects, e.g. `c("Date", "Tank")`. These
#'                       columns are converted to factors and stored in
#'                       metadata for the fitter to read.
#' @param duration_unit  Label for the unit of `duration`, stored in metadata.
#'                       Default `"hours"`.
#' @param temp_mean      Value to subtract from `temp` to form `temp_c`.
#'                       `NULL` (default) uses `mean(temp)`. Supply a fixed
#'                       value to align multiple datasets to a common centre.
#' @return A tibble with the standardised columns (plus `censored` and
#'         `survival_recorded` when `detection_limit` is given) and a
#'         `"tdt_meta"` attribute
#'         storing `temp_mean`, `duration_unit`, `random_effects`,
#'         `response_type` (`"count"` or `"proportion"`), `response_var`
#'         (the response column name for a proportion fit, else `NULL`),
#'         `proportion_eps` (the clamp half-width used for a proportion fit, else
#'         `NULL`), `detection_limit` (as supplied, else `NULL`), `censored`
#'         (`TRUE` when a `detection_limit` was supplied) and `n_censored` (the
#'         numbers of left- and right-censored rows, else `NULL`).
#' @examples
#' # Count data
#' raw <- data.frame(
#'   temperature_C = rep(c(30, 32, 34), each = 4),
#'   exposure_h    = rep(c(1, 2, 4, 8), times = 3),
#'   n             = 30L,
#'   alive         = c(29, 28, 25, 5, 30, 27, 18, 2, 28, 22, 10, 1)
#' )
#' standardize_data(raw,
#'                  temp     = "temperature_C",
#'                  duration = "exposure_h",
#'                  n_total  = "n",
#'                  n_surv   = "alive")
#'
#' # Continuous proportion (Beta) data
#' raw_p <- data.frame(
#'   temperature_C = rep(c(30, 32, 34), each = 4),
#'   exposure_h    = rep(c(1, 2, 4, 8), times = 3),
#'   fvfm_ratio    = c(0.95, 0.9, 0.7, 0.2, 0.92, 0.6, 0.3, 0, 0.8, 0.4, 0.1, 0)
#' )
#' standardize_data(raw_p,
#'                  temp       = "temperature_C",
#'                  duration   = "exposure_h",
#'                  proportion = "fvfm_ratio")
#'
#' # The instrument records 0 below 0.05: treat those zeros as censored
#' cens <- standardize_data(raw_p,
#'                          temp            = "temperature_C",
#'                          duration        = "exposure_h",
#'                          proportion      = "fvfm_ratio",
#'                          detection_limit = 0.05)
#' table(cens$censored)
#' @export
standardize_data <- function(data,
                             temp,
                             duration,
                             n_total        = NULL,
                             n_surv         = NULL,
                             n_dead         = NULL,
                             survival       = NULL,
                             mortality      = NULL,
                             proportion     = NULL,
                             proportion_eps = 0.001,
                             random_effects = NULL,
                             duration_unit  = "hours",
                             temp_mean      = NULL,
                             detection_limit = NULL) {

  response_args <- list(n_surv = n_surv, n_dead = n_dead, survival = survival,
                        mortality = mortality, proportion = proportion)
  if (sum(!vapply(response_args, is.null, logical(1))) != 1L) {
    stop("Supply exactly one of n_surv, n_dead, survival, mortality, or ",
         "proportion.", call. = FALSE)
  }
  is_proportion <- !is.null(proportion)
  if (!is_proportion && is.null(n_total)) {
    stop("Count responses (n_surv/n_dead/survival/mortality) require `n_total`.",
         call. = FALSE)
  }
  if (!is.null(detection_limit) && !is_proportion)
    stop("`detection_limit` applies to a continuous `proportion` response only ",
         "(a reading censored at an instrument limit); count responses are ",
         "exact.", call. = FALSE)
  if (is_proportion &&
      (!is.numeric(proportion_eps) || length(proportion_eps) != 1L ||
       !is.finite(proportion_eps) || proportion_eps <= 0 || proportion_eps >= 0.5))
    stop("`proportion_eps` must be one number in (0, 0.5).", call. = FALSE)

  needed <- c(temp, duration, n_total, n_surv, n_dead, survival, mortality,
              proportion, tdt_random_effect_variables(random_effects))
  tdt_check_columns(data, needed, "input columns")

  # Warn if standardising will clobber a pre-existing column whose name is a
  # transform-output name but which is NOT the source for that slot -- e.g. a
  # categorical column literally named `temp`, or a `temp_c`/`logd` that means
  # something else. (The count/response columns n_surv/n_dead/survival are
  # intentionally excluded: they are routinely present in count data and are
  # recomputed to consistent values, so warning on them is just noise.)
  reserved  <- c("temp", "duration", "logd", "temp_c",
                 if (!is.null(detection_limit)) c("censored", "survival_recorded"))
  sources   <- unlist(list(temp, duration, n_total, n_surv, n_dead, survival,
                           mortality, proportion,
                           Filter(is.character, as.list(detection_limit))))
  clobbered <- setdiff(intersect(reserved, names(data)), sources)
  if (length(clobbered))
    warning("standardize_data() overwrites existing column(s) with standardised ",
            "values: ", paste(clobbered, collapse = ", "),
            ". Rename them in the raw data to keep their original contents.",
            call. = FALSE)

  out          <- as.data.frame(data)
  out$temp     <- as.numeric(out[[temp]])
  out$duration <- as.numeric(out[[duration]])
  out$logd     <- log10(out$duration)

  if (is_proportion) {
    response_type <- "proportion"
    response_var  <- "survival"
    # The proportion branch stores its (clamped) values in `survival`. Unlike the
    # count paths -- where `survival` is a routine, recomputed response column --
    # a pre-existing `survival` here is unrelated user data that would be silently
    # destroyed, so warn unless it is itself the proportion source.
    if ("survival" %in% names(data) && !identical(proportion, "survival"))
      warning("standardize_data() overwrites the existing `survival` column with ",
              "the clamped proportion response. Rename it in the raw data to keep ",
              "its original contents.", call. = FALSE)
    p             <- as.numeric(out[[proportion]])
    if (any(p > 1 + 1e-6, na.rm = TRUE) || any(p < -1e-6, na.rm = TRUE))
      stop("`proportion` has values outside [0, 1]. A continuous-proportion ",
           "response must be a fraction in [0, 1] (e.g. post/pre Fv/Fm).",
           call. = FALSE)
    if (is.null(detection_limit)) {
      out$survival <- pmin(pmax(p, proportion_eps), 1 - proportion_eps)
    } else {
      # Censored readings: the limit is the response (brms reads it as the
      # bound of the censoring interval) and the flag says which side.
      # A reading AT the lower limit is a real reading (e.g. 0.081 / 0.81 is
      # 0.09999999999999999 in floating point, yet equals a limit of 0.1), and
      # one AT the upper limit is a capped reading, so compare with a small
      # relative tolerance on both sides.
      tol   <- sqrt(.Machine$double.eps)
      lim   <- tdt_resolve_detection_limit(detection_limit, out)
      left  <- !is.na(p) & !is.na(lim$lower) & p < lim$lower * (1 - tol)
      right <- !is.na(p) & !is.na(lim$upper) & p > lim$upper * (1 - tol)
      # Warn only about rows that survive the cleaning filter below.
      ok    <- is.finite(out$temp) & is.finite(out$duration) & out$duration > 0
      n_pos <- sum(left & p > lim$lower * tol & ok)     # a 0 (up to rounding) is expected
      if (n_pos > 0L)
        warning(n_pos, " positive reading(s) lie below the lower detection ",
                "limit and are treated as censored (true value below the ",
                "limit). If the instrument can read below this limit, lower it.",
                call. = FALSE)
      # A reading AT the upper limit (or of 1) is the expected capped value;
      # one strictly between the limit and 1 says the limit is set too low.
      n_sub1 <- sum(right & p > lim$upper * (1 + tol) & p < 1 - tol & ok)
      if (n_sub1 > 0L)
        warning(n_sub1, " reading(s) lie above the upper detection limit but ",
                "below 1 and are treated as censored (true value at or above ",
                "the limit). If the instrument can read above this limit, ",
                "raise it.", call. = FALSE)
      # Clamp only the sides without a limit: on a side with one, an
      # uncensored reading already lies strictly inside (0, 1).
      clamped <- pmin(pmax(p, ifelse(is.na(lim$lower), proportion_eps, 0)),
                      ifelse(is.na(lim$upper), 1 - proportion_eps, 1))
      out$survival_recorded <- p
      out$survival <- ifelse(left, lim$lower, ifelse(right, lim$upper, clamped))
      out$censored <- ifelse(left, "left", ifelse(right, "right", "none"))
      if (any(!is.na(p) & is.na(out$survival)))
        stop("Internal error: a reading lost its response value while applying ",
             "`detection_limit`. Please report this with the call.",
             call. = FALSE)
    }
    keep <- is.finite(out$survival) & is.finite(out$temp) &
            is.finite(out$duration) & out$duration > 0
  } else {
    response_type <- "count"
    response_var  <- NULL
    out$n_total   <- as.integer(round(as.numeric(out[[n_total]])))

    if (!is.null(n_surv)) {
      out$n_surv <- as.integer(round(as.numeric(out[[n_surv]])))
    } else if (!is.null(n_dead)) {
      out$n_surv <- out$n_total - as.integer(round(as.numeric(out[[n_dead]])))
    } else if (!is.null(survival)) {
      sv <- as.numeric(out[[survival]])
      if (any(sv > 1 + 1e-6, na.rm = TRUE))
        stop("`survival` has values > 1, which are not proportions. If these are ",
             "survivor COUNTS, pass them via `n_surv =` instead (with `n_total`).",
             call. = FALSE)
      out$n_surv <- as.integer(round(pmin(pmax(sv, 0), 1) * out$n_total))
    } else if (!is.null(mortality)) {
      mv <- as.numeric(out[[mortality]])
      if (any(mv > 1 + 1e-6, na.rm = TRUE))
        stop("`mortality` has values > 1, which are not proportions. If these are ",
             "death COUNTS, pass them via `n_dead =` instead (with `n_total`).",
             call. = FALSE)
      out$n_surv <- as.integer(round((1 - pmin(pmax(mv, 0), 1)) * out$n_total))
    }

    # Defensive clamp to [0, n_total]. The survival/mortality paths already
    # reject proportions > 1; the n_surv / n_dead count paths can still produce
    # out-of-range counts from data-entry errors (survivors > trials, deaths >
    # trials), so warn rather than silently fabricating a 0%/100% cell.
    n_oob <- sum(out$n_surv > out$n_total | out$n_surv < 0, na.rm = TRUE)
    if (n_oob > 0L)
      warning(n_oob, " cell(s) had survivor counts outside [0, n_total] ",
              "(e.g. n_surv > n_total); clamped to the valid range. ",
              "Check for data-entry errors.", call. = FALSE)
    out$n_surv   <- pmin(pmax(out$n_surv, 0), out$n_total)
    out$n_dead   <- out$n_total - out$n_surv
    out$survival <- out$n_surv / out$n_total
    keep <- is.finite(out$n_total) & is.finite(out$n_surv) &
            is.finite(out$temp)    & is.finite(out$duration) &
            out$n_total > 0 & out$duration > 0
  }

  out <- out[keep, , drop = FALSE]

  n_temp <- length(unique(out$temp))
  if (n_temp < 2L)
    warning("Only ", n_temp, " unique assay temperature(s) after cleaning: the ",
            "temperature slope (z = -1/slope) is not identified by the data and ",
            "would be driven entirely by the prior. Provide multiple assay ",
            "temperatures (>= 3 recommended) for an identified z.", call. = FALSE)

  if (is.null(temp_mean)) temp_mean <- mean(out$temp, na.rm = TRUE)
  out$temp_c <- out$temp - temp_mean

  for (re_var in tdt_random_effect_variables(random_effects)) {
    out[[re_var]] <- factor(out[[re_var]])
  }

  attr(out, "tdt_meta") <- list(
    temp_mean      = temp_mean,
    duration_unit  = duration_unit,
    random_effects = random_effects,
    response_type  = response_type,
    response_var   = response_var,
    proportion_eps = if (is_proportion) proportion_eps else NULL,
    detection_limit = detection_limit,
    censored       = !is.null(detection_limit),
    n_censored     = if (!is.null(detection_limit))
                       c(left  = sum(out$censored == "left"),
                         right = sum(out$censored == "right"))
  )

  tibble::as_tibble(out)
}

# Resolve `detection_limit` (see standardize_data()) into per-row lower and
# upper limits, NA where a side has no limit. Accepts a number (a lower limit),
# a length-2 numeric c(lower, upper) (named or positional; NA = no limit on that
# side), column name(s) holding per-row limits, or a list mixing the two.
# Every limit must lie strictly inside (0, 1) and below its upper partner.
tdt_resolve_detection_limit <- function(detection_limit, data) {
  x <- detection_limit
  if (!length(x) || length(x) > 2L)
    stop("`detection_limit` must have length 1 (a lower limit, or one named ",
         "`upper`) or 2 (c(lower, upper)).", call. = FALSE)
  nm <- names(x)
  if (is.null(nm) || all(nm == "")) {
    nm <- c("lower", "upper")[seq_along(x)]
  } else if (!all(nm %in% c("lower", "upper")) || anyDuplicated(nm)) {
    stop("Name the `detection_limit` entries `lower` and/or `upper`.",
         call. = FALSE)
  }
  n   <- nrow(data)
  lim <- list(lower = rep(NA_real_, n), upper = rep(NA_real_, n))
  for (i in seq_along(x)) {
    v <- x[[i]]
    if (length(v) != 1L)
      stop("Each `detection_limit` entry must be one number (or NA) or one ",
           "column name.", call. = FALSE)
    if (is.na(v)) {
      vals <- rep(NA_real_, n)              # NA of any type: no limit this side
    } else if (is.character(v)) {
      if (!nzchar(v) || !v %in% names(data))
        stop("`detection_limit` names a column that is not in the data: '", v,
             "'.", call. = FALSE)
      col <- data[[v]]
      if (!is.numeric(col) && !all(is.na(col)))
        stop("The `detection_limit` column `", v, "` must be numeric (it is ",
             class(col)[1], "); convert it first, checking entries such as ",
             "'<0.1' or '0,1'. Use NA for rows without a limit.", call. = FALSE)
      vals <- as.numeric(col)
    } else if (is.numeric(v)) {
      vals <- rep(as.numeric(v), n)
    } else {
      stop("Each `detection_limit` entry must be one number (or NA) or one ",
           "column name.", call. = FALSE)
    }
    stopifnot(length(vals) == n)
    if (any(!is.na(vals) & !(vals > 0 & vals < 1)))
      stop("The ", nm[i], " detection limit must lie strictly between 0 and 1 ",
           "(NA for no limit).", call. = FALSE)
    lim[[nm[i]]] <- vals
  }
  if (n > 0L && all(is.na(lim$lower)) && all(is.na(lim$upper)))
    stop("`detection_limit` sets no limit (every value is NA).", call. = FALSE)
  both <- !is.na(lim$lower) & !is.na(lim$upper)
  if (any(lim$lower[both] >= lim$upper[both]))
    stop("The lower detection limit must lie below the upper one.",
         call. = FALSE)
  lim
}
