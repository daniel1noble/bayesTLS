# Three-trace heat-injury validation harness:
#
#   - flat:         constant baseline, no damage accumulation expected.
#   - single_spike: one short pulse, analytical dose computable from z + CTmax.
#   - multi_spike:  several short pulses, dose = sum of individual contributions.
#
# Use planted_dose_from_trace() to compute the analytical truth given any trace
# and any (z, CTmax_1hr, T_c) triple.

#' Build reference temperature traces for heat-injury validation
#'
#' Returns a list of four tibbles, each with columns `time` and `temp`,
#' representing the canonical validation scenarios:
#'
#' - `flat` — constant `baseline` temperature for `n_hours` hours.
#' - `single_spike` — flat baseline with one time step (`dt_hours` wide) at
#'    `spike_temp` placed at each of `spike_times_single` (default one spike).
#' - `multi_spike` — flat baseline with one time step (`dt_hours` wide) at
#'    `spike_temp` placed at each of `spike_times_multi` (default three spikes).
#' - `diurnal` — multi-day diurnal cycle with day-to-day variation in peak
#'    temperature and small hourly fluctuations on top. Useful as a stand-in
#'    for a natural thermal regime where some days exceed `T_crit` and accrue
#'    injury while others stay below.
#'
#' Each spike is one time step (`dt_hours`) wide, so the analytical dose it
#' delivers is the instantaneous rate times the step width,
#' \eqn{100\cdot 10^{(T_{spike}-CT_{max,1hr})/z}\cdot dt_{hours}} % LT50-dose.
#' With the default `dt_hours = 1` this reduces to
#' \eqn{100\cdot 10^{(T_{spike}-CT_{max,1hr})/z}}, so setting
#' `spike_temp = CTmax_1hr` delivers ~100% LT50-dose per spike — a calibration
#' target for `predict_heat_injury()` validation. For `dt_hours != 1` the dose
#' scales by `dt_hours`.
#'
#' @param baseline    Numeric. Constant baseline temperature (°C). Default 20.
#' @param spike_temp  Numeric. Temperature of each spike (°C); each spike is one
#'                    time step (`dt_hours`) wide. Default 30.
#' @param n_hours     Integer. Total length of each trace, in hours. Default 96.
#' @param dt_hours    Numeric. Time step, in hours. Default 1 (hourly).
#' @param spike_times_single Integer vector. Hours at which the single-spike
#'                    trace has a spike. Default `24`.
#' @param spike_times_multi  Integer vector. Hours at which the multi-spike
#'                    trace has a spike. Default `c(24, 48, 72)`.
#' @param diurnal_n_days   Integer. Number of days in the diurnal scenario.
#'                    Default 7.
#' @param diurnal_night_temp Numeric vector of length `diurnal_n_days` (or
#'                    length 1, recycled): night-time baseline temperature (°C)
#'                    each day, applied at the cool end of that day's diurnal
#'                    cycle. Passing a vector lets nights on hot days be warmer
#'                    than nights on cool days (a common natural pattern).
#'                    Default `baseline` (constant across days).
#' @param diurnal_day_peaks  Numeric vector of length `diurnal_n_days` (or
#'                    length 1, recycled): the peak temperature each day reaches
#'                    near 14:00. If `NULL` (default) the function alternates
#'                    cool (`baseline + 5`) and warm (`baseline + 8`) days so
#'                    some days accrue HI and others do not.
#' @param diurnal_peak_fwhm Numeric. Full-width-half-maximum of the daily
#'                    Gaussian temperature peak, in hours. Default 6 (≈ 4-5
#'                    hours near the daily peak).
#' @param diurnal_noise_sd Numeric. Standard deviation (°C) of smooth hourly
#'                    fluctuations added on top of the diurnal cycle.
#'                    Default 0.3.
#' @param diurnal_seed Integer. Seed for the hourly-noise RNG so the diurnal
#'                    trace is reproducible. Default 1L.
#' @return A named list of four tibbles (`flat`, `single_spike`, `multi_spike`,
#'         `diurnal`), each with columns `time` (numeric, hours from start)
#'         and `temp` (°C). Note the traces differ in length: `flat`,
#'         `single_spike`, and `multi_spike` span `n_hours` (default 96) while
#'         `diurnal` spans `diurnal_n_days * 24` (default 168) hours.
#' @examples
#' scens <- make_temperature_scenarios()
#' lapply(scens, head)
#' # Diurnal calibrated for a shrimp-like organism (T_crit ~ 25 °C):
#' make_temperature_scenarios(
#'   baseline = 18,
#'   diurnal_n_days = 7,
#'   diurnal_night_temp = 19,
#'   diurnal_day_peaks  = c(24, 27, 24.5, 27.5, 23, 26.5, 25.5)
#' )$diurnal
#' @export
make_temperature_scenarios <- function(baseline    = 20,
                                       spike_temp  = 30,
                                       n_hours     = 96,
                                       dt_hours    = 1,
                                       spike_times_single = 24,
                                       spike_times_multi  = c(24, 48, 72),
                                       diurnal_n_days     = 7,
                                       diurnal_night_temp = baseline,
                                       diurnal_day_peaks  = NULL,
                                       diurnal_peak_fwhm  = 6,
                                       diurnal_noise_sd   = 0.3,
                                       diurnal_seed       = 1L) {

  time  <- seq(0, n_hours - dt_hours, by = dt_hours)
  base    <- tibble::tibble(time = time, temp = baseline)

  flat <- base

  apply_spikes <- function(spike_hours) {
    out <- base
    for (h in spike_hours) {
      i <- which(out$time == h)
      if (length(i) == 1L) out$temp[i] <- spike_temp
    }
    out
  }

  # ----- diurnal scenario -------------------------------------------------
  diurnal_t_h <- seq(0, diurnal_n_days * 24 - dt_hours, by = dt_hours)
  hour_of_day <- diurnal_t_h %% 24
  day_index   <- floor(diurnal_t_h / 24) + 1L

  if (is.null(diurnal_day_peaks)) {
    diurnal_day_peaks <- rep(c(baseline + 5, baseline + 8),
                              length.out = diurnal_n_days)
  } else if (length(diurnal_day_peaks) == 1L) {
    diurnal_day_peaks <- rep(diurnal_day_peaks, diurnal_n_days)
  } else if (length(diurnal_day_peaks) < diurnal_n_days) {
    diurnal_day_peaks <- rep(diurnal_day_peaks, length.out = diurnal_n_days)
  }

  # Per-day night baseline (vector). Recycle a length-1 scalar; pad shorter
  # vectors. Lets hot days couple to warmer nights.
  if (length(diurnal_night_temp) == 1L) {
    diurnal_night_temp <- rep(diurnal_night_temp, diurnal_n_days)
  } else if (length(diurnal_night_temp) < diurnal_n_days) {
    diurnal_night_temp <- rep(diurnal_night_temp, length.out = diurnal_n_days)
  }

  # Gaussian peak centred at 14:00 with the requested FWHM.
  sigma         <- diurnal_peak_fwhm / 2.355
  diurnal_shape <- exp(-((hour_of_day - 14)^2) / (2 * sigma^2))

  peak_at_t  <- diurnal_day_peaks[day_index]
  night_at_t <- diurnal_night_temp[day_index]
  diurnal_temp <- night_at_t +
                  (peak_at_t - night_at_t) * diurnal_shape

  # Smoothed hourly noise: random walk filtered with a length-3 moving average
  # so adjacent hours have correlated, "consistent" fluctuations rather than
  # white-noise jitter.
  set.seed(diurnal_seed)
  raw_noise   <- stats::rnorm(length(diurnal_t_h))
  smoothed    <- stats::filter(raw_noise, c(1, 2, 1) / 4, sides = 2)
  smoothed    <- as.numeric(smoothed)
  smoothed[is.na(smoothed)] <- 0
  if (stats::sd(smoothed) > 0) {
    smoothed <- smoothed * (diurnal_noise_sd / stats::sd(smoothed))
  }
  diurnal_temp <- diurnal_temp + smoothed

  diurnal <- tibble::tibble(time = diurnal_t_h, temp = diurnal_temp)

  list(
    flat         = flat,
    single_spike = apply_spikes(spike_times_single),
    multi_spike  = apply_spikes(spike_times_multi),
    diurnal      = diurnal
  )
}

#' Analytical heat-injury accumulation along a temperature trace
#'
#' Implements the classical HI integral exactly (Equation 8 of the manuscript),
#' given known TDT parameters. Use this as the **truth** that
#' [predict_heat_injury()] should recover when fed the same trace and a model
#' fit whose posterior is consistent with `(z, CTmax_1hr, T_c)`.
#'
#' Integrated by forward Euler from zero: the interval ending at point \eqn{i}
#' accrues the rate at its start (point \eqn{i-1}) times that interval's width,
#'
#' \deqn{\Delta HI_i = 100 \cdot 10^{(T_{i-1} - CT_{max,1hr}) / z} \cdot (t_i - t_{i-1})}
#'
#' when `T_{i-1} > T_c`, and zero otherwise (and `hi_inc[1] = 0`). Cumulative HI
#' is the running sum. This matches [predict_heat_injury()]'s integrator.
#'
#' @param trace      Tibble with columns `time` and `temp`, output of
#'                   [make_temperature_scenarios()].
#' @param z          Thermal sensitivity, °C per 10-fold change in time.
#' @param CTmax_1hr  Static temperature at which LT50 = 1 hour, °C.
#' @param T_c        Damage threshold (°C); contributions below `T_c` are zero.
#' @return The input trace augmented with two new columns:
#'         - `hi_inc` — instantaneous HI contribution at each hour (%).
#'         - `hi_cumulative` — running sum of `hi_inc` (%).
#' @examples
#' scens <- make_temperature_scenarios(spike_temp = 30)
#' planted_dose_from_trace(scens$single_spike, z = 5, CTmax_1hr = 32, T_c = 25)
#' @export
planted_dose_from_trace <- function(trace, z, CTmax_1hr, T_c) {
  if (anyNA(trace$time) || anyNA(trace$temp))
    stop("trace has NA in `time` or `temp`; interpolate or drop missing rows ",
         "before computing the planted dose.", call. = FALSE)
  n <- nrow(trace)
  # Instantaneous damage rate (%/hour) at each point; zero below the threshold.
  rate <- ifelse(trace$temp > T_c, 100 * 10 ^ ((trace$temp - CTmax_1hr) / z), 0)
  # Forward-Euler with per-interval widths, matching predict_heat_injury(): the
  # interval ending at point j accrues the rate at its START (point j-1) times
  # that interval's width; point 1 has accrued nothing yet. (The previous version
  # reused diff(time)[1] for every step and credited a step at t = 0, over-counting
  # by one step and mishandling irregular traces.)
  inc <- numeric(n)
  if (n >= 2L) inc[-1] <- rate[-n] * diff(trace$time)
  trace$hi_inc        <- inc
  trace$hi_cumulative <- cumsum(inc)
  trace
}
