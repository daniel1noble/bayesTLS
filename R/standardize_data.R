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
#'   finite); no `n_total`/`n_surv` columns are created. `response_type` is
#'   recorded as `"proportion"`.
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
#'                       at exactly 0 or 1). Default `0.001`.
#' @param random_effects Optional character vector of grouping variables for
#'                       random effects, e.g. `c("Date", "Tank")`. These
#'                       columns are converted to factors and stored in
#'                       metadata for the fitter to read.
#' @param duration_unit  Label for the unit of `duration`, stored in metadata.
#'                       Default `"hours"`.
#' @param temp_mean      Value to subtract from `temp` to form `temp_c`.
#'                       `NULL` (default) uses `mean(temp)`. Supply a fixed
#'                       value to align multiple datasets to a common centre.
#' @return A tibble with the standardised columns plus a `"tdt_meta"` attribute
#'         storing `temp_mean`, `duration_unit`, `random_effects`,
#'         `response_type` (`"count"` or `"proportion"`), and `response_var`
#'         (the response column name for a proportion fit, else `NULL`).
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
                             temp_mean      = NULL) {

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

  needed <- c(temp, duration, n_total, n_surv, n_dead, survival, mortality,
              proportion, tdt_random_effect_variables(random_effects))
  tdt_check_columns(data, needed, "input columns")

  out          <- as.data.frame(data)
  out$temp     <- as.numeric(out[[temp]])
  out$duration <- as.numeric(out[[duration]])
  out$logd     <- log10(out$duration)

  if (is_proportion) {
    response_type <- "proportion"
    response_var  <- "survival"
    p             <- as.numeric(out[[proportion]])
    out$survival  <- pmin(pmax(p, proportion_eps), 1 - proportion_eps)
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
      prop <- pmin(pmax(as.numeric(out[[survival]]), 0), 1)
      out$n_surv <- as.integer(round(prop * out$n_total))
    } else if (!is.null(mortality)) {
      prop <- 1 - pmin(pmax(as.numeric(out[[mortality]]), 0), 1)
      out$n_surv <- as.integer(round(prop * out$n_total))
    }

    out$n_surv   <- pmin(pmax(out$n_surv, 0), out$n_total)
    out$n_dead   <- out$n_total - out$n_surv
    out$survival <- out$n_surv / out$n_total
    keep <- is.finite(out$n_total) & is.finite(out$n_surv) &
            is.finite(out$temp)    & is.finite(out$duration) &
            out$n_total > 0 & out$duration > 0
  }

  out <- out[keep, , drop = FALSE]

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
    response_var   = response_var
  )

  tibble::as_tibble(out)
}
