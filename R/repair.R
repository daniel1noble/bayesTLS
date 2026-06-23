#' Sharpe-Schoolfield thermal performance curve for repair rate
#'
#' Computes a temperature-dependent repair rate using the Sharpe-Schoolfield
#' formulation (Sharpe & Schoolfield, 1981). All Arrhenius temperatures inside
#' this function are in Kelvin; the user-facing `temp_celsius` argument is in
#' degrees Celsius and converted internally.
#'
#' The Sharpe-Schoolfield rate at temperature `T_K` (in Kelvin) is
#'
#' \deqn{r(T_K) = \frac{r_{ref} \cdot \exp\bigl(T_A (T_{ref}^{-1} - T_K^{-1})\bigr)}
#'                    {1 + \exp\bigl(T_{AL} (T_K^{-1} - T_L^{-1})\bigr)
#'                       + \exp\bigl(T_{AH} (T_H^{-1} - T_K^{-1})\bigr)}}
#'
#' where the numerator is the Arrhenius enzymatic rate at the reference
#' temperature, and the two terms in the denominator suppress the rate at
#' temperatures below `TL` and above `TH` respectively (the low- and
#' high-temperature inactivation arms).
#'
#' @param temp_celsius Numeric vector of temperatures in °C.
#' @param TA  Arrhenius activation energy (K).
#' @param TAL Low-temperature inactivation activation energy (K).
#' @param TAH High-temperature inactivation activation energy (K).
#' @param TL  Low-temperature inactivation midpoint, in Kelvin.
#' @param TH  High-temperature inactivation midpoint, in Kelvin.
#' @param TREF Reference temperature, in Kelvin.
#' @param r_ref The Arrhenius (uninhibited) rate scale, in the same units the
#'              user wants the output expressed in. Note this is NOT exactly the
#'              realised rate at `TREF`: the inactivation denominator suppresses
#'              it, so `rate(TREF) = r_ref / (1 + exp(TAL(1/TREF - 1/TL)) +
#'              exp(TAH(1/TH - 1/TREF)))` (a few % below `r_ref` when `TREF` sits
#'              well between `TL` and `TH`, more as it approaches either). To set
#'              the realised rate at `TREF` to a target value, divide your target
#'              by that denominator.
#' @return Numeric vector of repair rates at the supplied temperatures, in the
#'         same units as `r_ref`. Negative or non-finite values are coerced to
#'         zero.
#' @examples
#' # Default shrimp-style parameters from the prototype: optimum ~17 °C
#' repair_rate_schoolfield(
#'   temp_celsius = seq(5, 30, by = 5),
#'   TA  = 14065, TAL = 50000, TAH = 120000,
#'   TL  = 10.5 + 273.15, TH = 22.5 + 273.15, TREF = 17 + 273.15,
#'   r_ref = 0.01
#' )
#' @export
repair_rate_schoolfield <- function(temp_celsius,
                                    TA, TAL, TAH,
                                    TL, TH, TREF,
                                    r_ref) {
  # TL/TH/TREF are Kelvin. Biological temperatures in Kelvin are ~250-330; a
  # value below 150 almost certainly means the caller passed Celsius (e.g.
  # TREF = 17 instead of 17 + 273.15), which silently collapses the rate to ~0
  # and turns repair off without warning. Catch it.
  k_args <- c(TL = TL, TH = TH, TREF = TREF)
  if (any(k_args < 150))
    stop("repair_rate_schoolfield(): TL/TH/TREF must be in KELVIN, but ",
         paste(names(k_args)[k_args < 150], collapse = ", "),
         " < 150 -- that looks like Celsius. Add 273.15 (e.g. TREF = 17 + 273.15).",
         call. = FALSE)
  temp_k <- temp_celsius + 273.15
  numerator   <- exp(TA * (1 / TREF - 1 / temp_k))
  denominator <- 1 +
    exp(TAL * (1 / temp_k - 1 / TL)) +
    exp(TAH * (1 / TH    - 1 / temp_k))
  rate <- numerator / denominator * r_ref
  rate[!is.finite(rate)] <- 0
  pmax(rate, 0)
}
