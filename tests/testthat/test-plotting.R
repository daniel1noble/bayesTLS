# Tests for the plotting helpers (R/plotting.R). Every plot function consumes
# a plain summary data frame (not a live brms fit), so these are brms-free: we
# build synthetic inputs whose shape matches the documented producer output,
# then assert via ggplot2::ggplot_build() that the *plotted geometry encodes
# the input data correctly* and that each option flips the expected scale /
# layer. This is stronger than "returns a ggplot": it checks the mapping.

# --- helpers ---------------------------------------------------------------
built <- function(p) ggplot2::ggplot_build(p)

# Does any scale for `aes` use a log-10 transform?
has_log_scale <- function(p, aes) {
  any(vapply(p$scales$scales, function(s) {
    aes_ok <- any(s$aesthetics %in% aes)
    tr <- s$trans %||% s$transform
    aes_ok && !is.null(tr) && identical(tr$name, "log-10")
  }, logical(1)))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ============================ theme_tdt ====================================

test_that("theme_tdt returns a theme honouring base_size and bold axis titles", {
  th <- theme_tdt(base_size = 20)
  expect_s3_class(th, "theme")
  expect_equal(th$text$size, 20)
  expect_equal(th$axis.title$face, "bold")
  # Attaches cleanly to a plot and builds.
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, hp)) +
    ggplot2::geom_point() + theme_tdt()
  expect_silent(built(p))
})

# ===================== plot_temperature_scenarios ==========================

test_that("plot_temperature_scenarios maps each scenario to its own facet and traces the temps", {
  scens <- make_temperature_scenarios(baseline = 20, spike_temp = 30,
                                      n_hours = 24, spike_times_single = 12,
                                      spike_times_multi = c(6, 12, 18))
  p  <- plot_temperature_scenarios(scens)
  bd <- built(p)
  # One facet panel per scenario.
  expect_equal(length(unique(bd$data[[1]]$PANEL)), length(scens))
  # The plotted y-values are exactly the concatenated scenario temperatures.
  plotted_y <- sort(bd$data[[1]]$y)
  raw_y     <- sort(unlist(lapply(scens, function(s) s$temp), use.names = FALSE))
  expect_equal(plotted_y, raw_y, tolerance = 1e-9)
})

test_that("plot_temperature_scenarios adds a T_c reference line only when asked", {
  scens <- make_temperature_scenarios(n_hours = 12)
  expect_equal(length(plot_temperature_scenarios(scens)$layers), 1L)
  p_tc <- plot_temperature_scenarios(scens, T_c = 24)
  expect_equal(length(p_tc$layers), 2L)            # line + hline
  # The hline sits at T_c.
  hline_layer <- p_tc$layers[[2]]
  expect_equal(hline_layer$data$yintercept, 24)
})

# ========================= plot_repair_rate ================================

test_that("plot_repair_rate plots exactly the Schoolfield repair rate", {
  rp <- list(TA = 14065, TAL = 50000, TAH = 120000,
             TL = 10.5 + 273.15, TH = 22.5 + 273.15,
             TREF = 17 + 273.15, r_ref = 0.01)
  grid  <- seq(5, 30, length.out = 50)
  p     <- plot_repair_rate(grid, rp)
  line_y <- built(p)$data[[1]]$y
  truth  <- repair_rate_schoolfield(grid, rp$TA, rp$TAL, rp$TAH,
                                    rp$TL, rp$TH, rp$TREF, rp$r_ref)
  expect_equal(line_y, truth, tolerance = 1e-9)    # the line IS the rate
})

# ====================== plot_temperature_density ===========================

# Synthetic extract_tdt()-style temperature posterior.
fake_temp_post <- function(center = 35, sd = 1.2, n = 4000) {
  set.seed(1)
  draws <- tibble::tibble(.draw = seq_len(n), temp = stats::rnorm(n, center, sd))
  q <- stats::quantile(draws$temp, c(0.025, 0.5, 0.975), names = FALSE)
  list(draws = draws,
       summary = tibble::tibble(temp_lower = q[1], temp_median = q[2],
                                temp_upper = q[3]))
}

test_that("plot_temperature_density draws the density and the CI segment at the summary bounds", {
  tp <- fake_temp_post()
  p  <- plot_temperature_density(tp)
  # Layers: density, CI segment, median point.
  expect_gte(length(p$layers), 3L)
  seg_layer <- p$layers[[2]]
  # The credible-interval segment spans temp_lower -> temp_upper.
  expect_equal(seg_layer$data$temp_lower, tp$summary$temp_lower)
  expect_equal(seg_layer$data$temp_upper, tp$summary$temp_upper)
  expect_silent(built(p))
})

test_that("plot_temperature_density adds a truth line only when supplied", {
  tp <- fake_temp_post()
  expect_equal(length(plot_temperature_density(tp)$layers), 3L)
  p_truth <- plot_temperature_density(tp, truth = 36)
  expect_equal(length(p_truth$layers), 4L)
  expect_equal(p_truth$layers[[4]]$data$xintercept %||%
               p_truth$layers[[4]]$geom_params$xintercept, 36)
})

# ========================= plot_heat_injury ================================

fake_hi <- function() {
  t <- 0:48
  hi <- pmin(150, 3 * pmax(0, t - 6))
  list(summary = tibble::tibble(
    time = t,
    hi_median = hi, hi_lower = hi * 0.8, hi_upper = hi * 1.2,
    surv_median = pmax(0, 1 - hi / 200),
    surv_lower  = pmax(0, 1 - hi / 160),
    surv_upper  = pmin(1, 1 - hi / 260)))
}

test_that("plot_heat_injury returns two stacked panels with the HI threshold line", {
  hi <- fake_hi()
  p  <- plot_heat_injury(hi, lt50_threshold = 100)
  expect_s3_class(p, "patchwork")
  # Two sub-plots: HI (top) and survival (bottom).
  expect_equal(length(p$patches$plots) + 1L, 2L)   # patchwork holds n-1 + base
  # Rebuild the HI sub-plot directly to confirm the threshold + data.
  p_hi <- plot_heat_injury(hi)[[1]]
  expect_equal(built(p_hi)$data[[2]]$y, hi$summary$hi_median, tolerance = 1e-9)
})

# ====================== plot_survival_curves ===============================

fake_pred <- function() {
  grid <- expand.grid(temp = c(30, 34, 38), duration = c(0.5, 1, 2, 4))
  grid$survival_median <- plogis(2 - 0.5 * grid$duration - 0.1 * (grid$temp - 34))
  grid$survival_lower  <- pmax(0, grid$survival_median - 0.1)
  grid$survival_upper  <- pmin(1, grid$survival_median + 0.1)
  list(summary = tibble::as_tibble(grid))
}

test_that("plot_survival_curves encodes the survival summary and one line per temperature", {
  pr <- fake_pred()
  p  <- plot_survival_curves(pr)
  bd <- built(p)
  # Ribbon (layer 1) ymin/ymax track the credible band.
  expect_equal(sort(bd$data[[1]]$ymin), sort(pr$summary$survival_lower),
               tolerance = 1e-9)
  # One colour group per temperature.
  expect_equal(length(unique(bd$data[[2]]$colour)),
               length(unique(pr$summary$temp)))
})

test_that("plot_survival_curves log_time flips the duration axis to log10", {
  pr <- fake_pred()
  expect_false(has_log_scale(plot_survival_curves(pr, log_time = FALSE), "x"))
  expect_true(has_log_scale(plot_survival_curves(pr, log_time = TRUE), "x"))
})

test_that("plot_survival_curves overlays observed points only when provided", {
  pr  <- fake_pred()
  obs <- tibble::tibble(temp = c(30, 34), duration = c(1, 2),
                        survival = c(0.8, 0.4))
  expect_equal(length(plot_survival_curves(pr)$layers), 2L)        # ribbon+line
  expect_equal(length(plot_survival_curves(pr, observed = obs)$layers), 3L)
})

# ========================= plot_tdt_landscape ==============================

fake_landscape <- function() {
  # Evenly spaced grid (geom_raster wants even intervals) with survival that
  # spans the full 0-1 range, so the 0.25/0.5/0.75 contours all exist.
  grid <- expand.grid(temp = seq(30, 38, length.out = 12),
                      duration = seq(0.5, 30, length.out = 12))
  grid$survival_median <- plogis(8 - 0.4 * grid$duration - 0.5 * (grid$temp - 30))
  list(summary = tibble::as_tibble(grid))
}

test_that("plot_tdt_landscape builds a raster heatmap of survival with contours", {
  lsp <- fake_landscape()
  p   <- plot_tdt_landscape(lsp)
  expect_s3_class(p$layers[[1]]$geom, "GeomRaster")
  expect_s3_class(p$layers[[2]]$geom, "GeomContour")
  # geom_contour over a fill raster emits a benign "fill dropped" message (the
  # contour stat ignores fill); it occurs in real use too and is not a bug.
  bd <- suppressWarnings(built(p))
  # One raster cell per grid row.
  expect_equal(nrow(bd$data[[1]]), nrow(lsp$summary))
  # Contours are drawn at exactly the requested survival levels (and they exist
  # because the synthetic surface spans 0-1).
  expect_setequal(unique(bd$data[[2]]$level), c(0.25, 0.5, 0.75))
  # Survival fill scale is clamped to the [0, 1] probability range.
  fill_scale <- p$scales$scales[[which(vapply(p$scales$scales,
                  function(s) any(s$aesthetics == "fill"), logical(1)))[1]]]
  expect_equal(fill_scale$limits, c(0, 1))
})

test_that("plot_tdt_landscape log_time flips the duration axis and overlays observed", {
  lsp <- fake_landscape()
  expect_false(has_log_scale(plot_tdt_landscape(lsp), "y"))
  expect_true(has_log_scale(plot_tdt_landscape(lsp, log_time = TRUE), "y"))
  obs <- tibble::tibble(temp = 34, duration = 1, survival = 0.5)
  expect_gt(length(plot_tdt_landscape(lsp, observed = obs)$layers),
            length(plot_tdt_landscape(lsp)$layers))
})

# ========================== plot_tdt_curve =================================

fake_ltx <- function(target = "p=0.500", unit = "min", big = FALSE) {
  temp <- seq(30, 38, by = 1)
  med  <- 10^(3 - 0.15 * (temp - 30))            # falls with temperature
  upr  <- med * if (big) 50 else 1.5             # `big` pushes past the 48 h cap
  list(summary = tibble::tibble(
         target_surv = target, temp = temp,
         duration_lower = med * 0.7, duration_median = med,
         duration_upper = upr),
       output_time_unit = unit)
}

test_that("plot_tdt_curve returns two panels by default and single panels on request", {
  ltx <- fake_ltx()
  expect_s3_class(plot_tdt_curve(ltx), "patchwork")              # both
  expect_s3_class(plot_tdt_curve(ltx, panels = "linear"), "ggplot")
  p_log <- plot_tdt_curve(ltx, panels = "log")
  expect_s3_class(p_log, "ggplot")
  expect_true(has_log_scale(p_log, "y"))                        # log panel
})

test_that("plot_tdt_curve renders the target-survival label for each label style", {
  expect_match(plot_tdt_curve(fake_ltx("p=0.500"), panels = "linear")$labels$y,
               "50% survival")
  expect_match(plot_tdt_curve(fake_ltx("(low+up)/2"), panels = "linear")$labels$y,
               "low \\+ up")
  expect_match(plot_tdt_curve(fake_ltx("0.1"), panels = "linear")$labels$y,
               "10% survival")
})

test_that("plot_tdt_curve clamps the y-axis to 48 h only when the band exceeds it", {
  # Minutes: 48 h cap = 2880. `big = TRUE` pushes duration_upper well past it.
  capped <- plot_tdt_curve(fake_ltx(unit = "min", big = TRUE), panels = "linear")
  expect_equal(capped$coordinates$limits$y, c(0, 48 * 60))
  # Band stays under the cap -> no coord clamp.
  free <- plot_tdt_curve(fake_ltx(unit = "min", big = FALSE), panels = "linear")
  expect_null(free$coordinates$limits$y)
})

test_that("plot_tdt_curve carries the output time unit into the axis label", {
  expect_match(plot_tdt_curve(fake_ltx(unit = "hours"), panels = "linear")$labels$y,
               "hours")
})
