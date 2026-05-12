#' Standardise a raw survival / proportion dataset for the TDT function library
#'
#' Rewrites user column names into a single project-standard schema
#' (`temp`, `duration`, `logd`, `temp_c`, `n_total`, `n_surv`, `n_dead`,
#' `survival`) and attaches metadata used by every downstream fitting and
#' prediction helper. This is the single entry point for raw data — everything
#' else in the library assumes the output of this function.
#'
#' Supply **exactly one** of `n_surv`, `n_dead`, `survival`, or `mortality`.
#' The other counts are derived.
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
#' @param n_surv         Column name for survivor counts.
#' @param n_dead         Column name for death counts. Converted to `n_surv`
#'                       via `n_surv = n_total - n_dead`.
#' @param survival       Column name for survival proportions in `[0, 1]`.
#'                       Converted to integer counts via `n_total`.
#' @param mortality      Column name for mortality proportions in `[0, 1]`.
#'                       Converted to `n_surv = round((1 - mortality) * n_total)`.
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
#'         storing `temp_mean`, `duration_unit`, and `random_effects`.
#' @examples
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
#' @export
standardize_data <- function(data,
                             temp,
                             duration,
                             n_total,
                             n_surv         = NULL,
                             n_dead         = NULL,
                             survival       = NULL,
                             mortality      = NULL,
                             random_effects = NULL,
                             duration_unit  = "hours",
                             temp_mean      = NULL) {

  needed <- c(temp, duration, n_total, n_surv, n_dead, survival, mortality,
              tdt_random_effect_variables(random_effects))
  tdt_check_columns(data, needed, "input columns")

  count_args <- list(n_surv = n_surv, n_dead = n_dead,
                     survival = survival, mortality = mortality)
  if (sum(!vapply(count_args, is.null, logical(1))) != 1L) {
    stop("Supply exactly one of n_surv, n_dead, survival, or mortality.",
         call. = FALSE)
  }

  out <- as.data.frame(data)
  out$n_total <- as.integer(round(as.numeric(out[[n_total]])))

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
  out$temp     <- as.numeric(out[[temp]])
  out$duration <- as.numeric(out[[duration]])
  out$logd     <- log10(out$duration)

  keep <- is.finite(out$n_total) & is.finite(out$n_surv) &
          is.finite(out$temp)    & is.finite(out$duration) &
          out$n_total > 0 & out$duration > 0
  out <- out[keep, , drop = FALSE]

  if (is.null(temp_mean)) temp_mean <- mean(out$temp, na.rm = TRUE)
  out$temp_c <- out$temp - temp_mean

  for (re_var in tdt_random_effect_variables(random_effects)) {
    out[[re_var]] <- factor(out[[re_var]])
  }

  attr(out, "tdt_meta") <- list(
    temp_mean      = temp_mean,
    duration_unit  = duration_unit,
    random_effects = random_effects
  )

  tibble::as_tibble(out)
}
