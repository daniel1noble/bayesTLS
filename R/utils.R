# Small utilities used across the TDT function library.
# Internal helpers carry @keywords internal.

#' Quantile wrapper with TDT-friendly defaults
#'
#' @param x Numeric vector.
#' @param probs Numeric vector of quantile probabilities.
#' @return Numeric vector of length `length(probs)`.
#' @examples
#' tdt_quantile(rnorm(100))
#' @export
tdt_quantile <- function(x, probs = c(0.025, 0.5, 0.975)) {
  stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
}

#' Format a posterior median plus credible interval as a single string
#'
#' @param median,lower,upper Numeric scalars.
#' @param digits Integer rounding precision.
#' @return Character scalar like `"5.12 [4.87, 5.40]"`.
#' @examples
#' format_interval(5.123, 4.872, 5.401)
#' @export
format_interval <- function(median, lower, upper, digits = 2) {
  paste0(round(median, digits), " [",
         round(lower, digits), ", ",
         round(upper, digits), "]")
}

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Convert various clock formats to minutes
#'
#' Accepts POSIXt, hms / difftime, numeric fractions of a day (Excel time),
#' bare numerics (assumed minutes), or `"HH:MM:SS"` character strings.
#'
#' @param x Time value(s).
#' @return Numeric vector of minutes since 00:00.
#' @examples
#' clock_to_minutes("08:30:00")
#' clock_to_minutes(0.5)     # half a day = 720 min
#' @export
clock_to_minutes <- function(x) {
  if (inherits(x, "POSIXt")) {
    return(as.numeric(format(x, "%H")) * 60 +
             as.numeric(format(x, "%M")) +
             as.numeric(format(x, "%S")) / 60)
  }
  if (inherits(x, "hms") || inherits(x, "difftime")) {
    return(as.numeric(x, units = "mins"))
  }
  if (is.numeric(x)) {
    # Excel clock times often arrive as fractions of a day in [0, 1].
    if (all(is.na(x) | (x >= 0 & x <= 1))) return(x * 24 * 60)
    return(x)
  }
  parsed <- suppressWarnings(as.POSIXct(as.character(x),
                                        format = "%H:%M:%S", tz = "UTC"))
  as.numeric(format(parsed, "%H")) * 60 +
    as.numeric(format(parsed, "%M")) +
    as.numeric(format(parsed, "%S")) / 60
}

# --- internal helpers ---------------------------------------------------------

#' Error on missing columns
#'
#' @keywords internal
tdt_check_columns <- function(data, cols, arg_name = "columns") {
  cols <- cols[!is.na(cols) & nzchar(cols)]
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop("Missing ", arg_name, ": ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Format bare names as `(1 | name)` random-effect terms
#'
#' @keywords internal
tdt_format_random_effects <- function(random_effects = NULL) {
  if (is.null(random_effects) || length(random_effects) == 0) return(character())
  out <- vapply(random_effects, function(term) {
    term <- trimws(term)
    if (grepl("^\\(", term)) term else paste0("(1 | ", term, ")")
  }, character(1))
  unname(out)
}

#' Extract variable names from random-effect terms
#'
#' @keywords internal
tdt_random_effect_variables <- function(random_effects = NULL) {
  terms <- tdt_format_random_effects(random_effects)
  unique(unlist(lapply(terms, function(term) {
    vars <- all.vars(stats::as.formula(paste("~", term)))
    vars[vars != "1"]
  })))
}

#' Derive 4PL asymptote intervals from a user-supplied response range
#'
#' Given the lower and upper bounds of where the asymptotes can sit, returns
#' the disjoint intervals used by [make_4pl_formula()]'s `inv_logit` reparam.
#' `low` is mapped to `(lower + pad, midpoint - gap/2)`, `up` to
#' `(midpoint + gap/2, upper - pad)`. The gap kills label-switching by ensuring
#' `up > low` always; the pad keeps the asymptotes off the exact boundaries.
#'
#' @param lower,upper Numeric scalars. The response-scale range that the
#'                    asymptotes can occupy (`0` and `1` for proportion data;
#'                    `0.85` and `1` for PSII-like sublethal data, etc.).
#' @param pad Absolute padding from `lower` and `upper`. Default `0.001`.
#' @param gap Absolute gap between the low and up intervals. Default `0.002`.
#' @return Named list with `low_min`, `low_max`, `low_w`, `up_min`, `up_max`,
#'         `up_w`, `midpoint`.
#' @keywords internal
compute_4pl_bounds <- function(lower = 0, upper = 1,
                               pad = 0.001, gap = 0.002) {
  if (upper <= lower)
    stop("upper must be strictly greater than lower.", call. = FALSE)
  if (2 * pad + gap >= (upper - lower))
    stop("pad and gap leave no room for asymptote intervals; ",
         "reduce pad/gap or widen lower/upper.", call. = FALSE)

  midpoint <- (lower + upper) / 2
  low_min  <- lower + pad
  low_max  <- midpoint - gap / 2
  up_min   <- midpoint + gap / 2
  up_max   <- upper - pad

  list(low_min  = low_min,
       low_max  = low_max,
       low_w    = low_max - low_min,
       up_min   = up_min,
       up_max   = up_max,
       up_w     = up_max - up_min,
       midpoint = midpoint)
}
