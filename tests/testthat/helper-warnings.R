# Muffle only tls()'s / derive_temperature_for_duration()'s "no single
# crossing" warning (draws whose bent LT curve has no single crossing of its
# target within the temperatures searched), so a test that compares finite
# draws still surfaces every other warning.
muffle_uncrossed <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("no single crossing", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
}
