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
#' @param median,lower,upper Numeric (scalar or vector).
#' @param digits Integer rounding precision.
#' @return Character like `"5.12 [4.87, 5.4]"`. A non-finite `median` yields
#'   `NA_character_` (rather than `"NA [...]"`); a non-finite bound is shown as an
#'   en dash, so the strings stay table-ready.
#' @examples
#' format_interval(5.123, 4.872, 5.401)
#' format_interval(NA, 1, 2)   # -> NA
#' @export
format_interval <- function(median, lower, upper, digits = 2) {
  fmt <- function(x) ifelse(is.finite(x), as.character(round(x, digits)), "\u2013")
  out <- paste0(fmt(median), " [", fmt(lower), ", ", fmt(upper), "]")
  out[!is.finite(median)] <- NA_character_
  out
}

#' Null-coalescing operator
#'
#' Returns `x` unless it is `NULL`, in which case `y` is returned.
#'
#' @param x,y Any R objects.
#' @return `x` if non-`NULL`, otherwise `y`.
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Convert various clock formats to minutes
#'
#' Accepts POSIXt, hms / difftime, numeric fractions of a day (Excel time) or
#' bare numeric minutes, and character strings: `"HH:MM:SS"`, `"HH:MM"`, bare
#' numeric strings (minutes), and durations beyond 24 h (e.g. `"25:30:00"`).
#' Character strings are parsed element-wise; malformed entries become `NA`.
#'
#' @param x Time value(s).
#' @return Numeric vector of minutes.
#' @examples
#' clock_to_minutes("08:30:00")
#' clock_to_minutes("25:30")   # 25 h 30 min = 1530 min
#' clock_to_minutes(0.5)       # half a day = 720 min
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
    nona <- x[!is.na(x)]
    # Excel clock times arrive as fractions of a day in [0, 1]; everything else
    # is taken as already-minutes. A vector that mixes the two is ambiguous --
    # warn instead of silently reinterpreting it (the old all-or-nothing rule
    # turned c(0.5, 720) into 0.5 min, not 720 min).
    if (length(nona) && all(nona >= 0 & nona <= 1)) {
      message("clock_to_minutes(): all numeric values are in [0, 1]; ",
              "treating them as Excel day-fractions (\u00d7 1440 min). ",
              "If they are already minutes, pass them as values > 1 or ",
              "convert with as.numeric() first.")
      return(x * 24 * 60)
    }
    if (length(nona) && any(nona > 0 & nona < 1) && any(nona > 1))
      warning("clock_to_minutes(): numeric input mixes values < 1 and > 1; ",
              "treating ALL as minutes. If the < 1 values are Excel day-fractions ",
              "they will be wrong -- convert the column to one unit first.",
              call. = FALSE)
    return(x)
  }
  # Character: parse element-wise so HH:MM, HH:MM:SS, durations > 24 h, and bare
  # numeric strings all work (as.POSIXct("%H:%M:%S") returned NA for all three).
  vapply(as.character(x), function(s) {
    if (is.na(s)) return(NA_real_)
    s <- trimws(s)
    if (!nzchar(s)) return(NA_real_)
    if (grepl("^[0-9]*\\.?[0-9]+$", s)) return(as.numeric(s))      # bare -> minutes
    parts <- strsplit(s, ":", fixed = TRUE)[[1]]
    if (!length(parts) %in% c(2L, 3L) ||
        any(!grepl("^[0-9]*\\.?[0-9]+$", parts))) return(NA_real_)
    h <- as.numeric(parts[1]); m <- as.numeric(parts[2])
    sec <- if (length(parts) == 3L) as.numeric(parts[3]) else 0
    h * 60 + m + sec / 60                                          # supports h >= 24
  }, numeric(1), USE.NAMES = FALSE)
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

#' Convert a time-unit label to minutes
#'
#' Maps a free-text duration/time unit (e.g. `"hours"`, `"min"`, `"s"`) to its
#' length in minutes. Used to derive the model-to-output `time_multiplier` in
#' [extract_tdt()] and [derive_tdt_curve()] from a workflow's `duration_unit`.
#'
#' @param unit Character scalar time-unit label.
#' @return Numeric scalar: the unit's length in minutes.
#' @keywords internal
tdt_unit_to_minutes <- function(unit) {
  if (is.null(unit) || length(unit) != 1L || is.na(unit))
    stop("time unit is NULL/NA", call. = FALSE)
  u <- tolower(trimws(as.character(unit)))
  switch(u,
         "s" = , "sec" = , "secs" = , "second" = , "seconds" = 1 / 60,
         "m" = , "min" = , "mins" = , "minute" = , "minutes" = 1,
         "h" = , "hr" = , "hrs" = , "hour" = , "hours"   = 60,
         "d" = , "day" = , "days" = 1440,
         stop("Unrecognised time unit: '", unit, "'. Use seconds/minutes/",
              "hours/days (or pass time_multiplier explicitly).", call. = FALSE))
}

#' Resolve the model-to-output time multiplier for TDT helpers
#'
#' If `time_multiplier` is supplied it is returned unchanged (explicit override).
#' Otherwise it is derived from the workflow's `meta$duration_unit` (the unit of
#' the model's `duration` column) and the requested `output_time_unit`, so that
#' `model_time * time_multiplier` is in `output_time_unit`. Falls back to `1`
#' (with a message) when the units cannot be resolved.
#'
#' @param time_multiplier Numeric scalar or `NULL`.
#' @param meta A `bayes_tls` workflow's `meta` list (uses `duration_unit`).
#' @param output_time_unit Target output time unit (e.g. `"min"`).
#' @return Numeric scalar multiplier.
#' @keywords internal
tdt_resolve_time_multiplier <- function(time_multiplier, meta,
                                        output_time_unit) {
  if (!is.null(time_multiplier)) return(time_multiplier)

  du <- meta$duration_unit
  mins_model  <- tryCatch(tdt_unit_to_minutes(du),
                          error = function(e) NA_real_)
  mins_output <- tryCatch(tdt_unit_to_minutes(output_time_unit),
                          error = function(e) NA_real_)

  if (is.na(mins_model) || is.na(mins_output)) {
    message("Could not derive time_multiplier from duration_unit ('",
            du %||% "NULL", "') and output_time_unit ('", output_time_unit,
            "'); assuming model time unit == output unit (time_multiplier = 1). ",
            "Pass time_multiplier explicitly to override.")
    return(1)
  }
  mins_model / mins_output
}

# --- grouped-fit detection ---------------------------------------------------
# A direct/midpoint fit is "grouped" when CTmax/z (or mid) vary by a fixed
# moderator. fit_4pl() records this coding-independently in meta$grouped, so
# `~ 0 + G` (cell-means) and `~ 1 + G` (treatment) are both detected. The
# coefficient-count fallback covers hand-wrapped brmsfits lacking that metadata.
# Used by tdt_resolve_by() as a safety net: a fit that looks grouped but has no
# recorded moderators must not be silently extracted at its reference level.
tdt_is_grouped <- function(workflow) {
  g <- workflow$meta$grouped
  if (!is.null(g)) return(isTRUE(g))
  if (!inherits(workflow, "bayes_tls") || is.null(workflow$fit)) return(FALSE)
  nms <- names(posterior::as_draws_df(workflow$fit))
  ndirect <- max(sum(grepl("^b_CTmaxdev_", nms)), sum(grepl("^b_logz_", nms)))
  if (ndirect > 0) return(ndirect > 1)
  any(grepl("^b_mid_", nms) & !grepl("_(Intercept|temp_c)$", nms))
}
