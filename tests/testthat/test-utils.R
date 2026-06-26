test_that("tdt_quantile returns the right length and values", {
  x <- rnorm(1000)
  q <- tdt_quantile(x)
  expect_length(q, 3)
  expect_equal(q[2], median(x), tolerance = 0.05)
})

test_that("format_interval rounds and formats correctly", {
  expect_equal(format_interval(5.123, 4.872, 5.401),
               "5.12 [4.87, 5.4]")
  expect_equal(format_interval(5.123, 4.872, 5.401, digits = 1),
               "5.1 [4.9, 5.4]")
})

test_that("%||% returns y when x is NULL", {
  expect_equal(NULL %||% 1, 1)
  expect_equal(2 %||% 1, 2)
})

test_that("compute_4pl_bounds gives the expected proportion defaults", {
  b <- compute_4pl_bounds(0, 1)
  expect_equal(b$low_min,  0.001)
  expect_equal(b$low_max,  0.499)
  expect_equal(b$up_min,   0.501)
  expect_equal(b$up_max,   0.999)
  expect_equal(b$midpoint, 0.5)
})

test_that("compute_4pl_bounds scales correctly for PSII-like ranges", {
  b <- compute_4pl_bounds(0.85, 1)
  expect_equal(b$midpoint, 0.925)
  expect_true(b$low_max < b$up_min)
  expect_true(b$up_max < 1)
  expect_true(b$low_min > 0.85)
})

test_that("compute_4pl_bounds errors when pad/gap exceed the range", {
  expect_error(compute_4pl_bounds(0.99, 1, pad = 0.01),
               "no room for asymptote intervals")
  expect_error(compute_4pl_bounds(1, 0),
               "upper must be strictly greater than lower")
})

test_that("tdt_format_random_effects wraps bare names but preserves explicit terms", {
  expect_equal(tdt_format_random_effects(c("Date", "Tank")),
               c("(1 | Date)", "(1 | Tank)"))
  expect_equal(tdt_format_random_effects("(1 + temp_c | Tank)"),
               "(1 + temp_c | Tank)")
  expect_equal(tdt_format_random_effects(NULL), character())
})

test_that("tdt_unit_to_minutes maps common time units", {
  expect_equal(tdt_unit_to_minutes("minutes"), 1)
  expect_equal(tdt_unit_to_minutes("min"),     1)
  expect_equal(tdt_unit_to_minutes("hours"),   60)
  expect_equal(tdt_unit_to_minutes("h"),       60)
  expect_equal(tdt_unit_to_minutes("seconds"), 1 / 60)
  expect_equal(tdt_unit_to_minutes("days"),    1440)
  expect_error(tdt_unit_to_minutes("fortnights"), "Unrecognised time unit")
})

test_that("tdt_resolve_time_multiplier derives from duration_unit, override honoured", {
  # Minutes model -> output min: multiplier 1 (the leaf/seed case the bug hit).
  expect_equal(
    tdt_resolve_time_multiplier(NULL, list(duration_unit = "minutes"), "min"), 1)
  # Hours model -> output min: multiplier 60 (shrimp/zebrafish).
  expect_equal(
    tdt_resolve_time_multiplier(NULL, list(duration_unit = "hours"), "min"), 60)
  # Explicit value overrides the derivation.
  expect_equal(
    tdt_resolve_time_multiplier(5, list(duration_unit = "hours"), "min"), 5)
  # Unknown unit falls back to 1 with a message.
  expect_message(
    out <- tdt_resolve_time_multiplier(NULL, list(duration_unit = NULL), "min"),
    "Could not derive time_multiplier")
  expect_equal(out, 1)
})

test_that("tdt_resolve_t_ref inherits the fit's reference (with a message) unless supplied", {
  # Explicit t_ref is returned unchanged and silently (the common case).
  expect_silent(out <- tdt_resolve_t_ref(60, list(t_ref = 240)))
  expect_equal(out, 60)
  # Omitted t_ref inherits the fit's stored reference and announces it -- this is
  # the footgun that read a 4 h fit at the 60 min default (CTmax off by z*log10(4)).
  expect_message(out <- tdt_resolve_t_ref(NULL, list(t_ref = 240)),
                 "model fit's reference exposure: 240 min \\(4 hours\\)")
  expect_equal(out, 240)
  # Raw fit with no recorded reference falls back to 60, also with a message.
  expect_message(out <- tdt_resolve_t_ref(NULL, list()), "defaulting to 60 min")
  expect_equal(out, 60)
})

test_that("tdt_format_hours renders singular/plural and fractional hours", {
  expect_equal(tdt_format_hours(60),  "1 hour")
  expect_equal(tdt_format_hours(240), "4 hours")
  expect_equal(tdt_format_hours(90),  "1.5 hours")
})

test_that("clock_to_minutes handles common formats", {
  expect_equal(clock_to_minutes("08:30:00"), 510)
  expect_equal(suppressMessages(clock_to_minutes(0.5)), 720)
  expect_equal(clock_to_minutes(as.POSIXct("2026-01-01 08:30:00", tz = "UTC")),
               510)
})

test_that("clock_to_minutes parses HH:MM, >24h, bare numbers, and NA element-wise", {
  expect_equal(clock_to_minutes("08:30:00"), 510)
  expect_equal(clock_to_minutes("08:30"),    510)   # HH:MM (no seconds)
  expect_equal(clock_to_minutes("25:30:00"), 1530)  # > 24 h
  expect_equal(clock_to_minutes("90"),       90)    # bare number -> minutes
  expect_equal(suppressMessages(clock_to_minutes(0.5)), 720) # Excel day-fraction
  expect_equal(clock_to_minutes(c("01:00:00", "02:00:00", "bad")),
               c(60, 120, NA))                       # malformed -> NA, not all-NA
})

test_that("clock_to_minutes warns (not silently flips) on an ambiguous mixed numeric vector", {
  expect_warning(r <- clock_to_minutes(c(0.5, 720)), "mixes values")
  expect_equal(r, c(0.5, 720))   # treated as minutes, the < 1 value left as-is
})

test_that("clock_to_minutes messages when applying the all-in-[0,1] day-fraction rule", {
  expect_message(r <- clock_to_minutes(c(0.25, 0.5)), "day-fractions")
  expect_equal(r, c(360, 720))
  # > 1 values are unambiguous minutes and stay silent.
  expect_silent(clock_to_minutes(c(60, 120)))
})
