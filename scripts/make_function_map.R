#!/usr/bin/env Rscript

# Draw a compact map of the bayesTLS user-facing workflow.
# The map intentionally omits low-level engine helpers and one-line wrappers.

out_dir <- file.path("output", "figs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cols <- list(
  ink = "#24313A",
  muted = "#64727D",
  grid = "#D7E0E5",
  data = "#DDEFE9",
  fit = "#DEE9F6",
  extract = "#E3F1DA",
  predict = "#F8E8CF",
  output = "#EEE5F4",
  plot = "#F2DEE8",
  optional = "#F1F3F4",
  warn = "#F5E0D8"
)

wrap <- function(x, width = 23) {
  paste(strwrap(x, width = width), collapse = "\n")
}

box <- function(x, y, w, h, title, body = NULL, fill = "white",
                border = cols$ink, title_cex = 0.78, body_cex = 0.58,
                lty = 1, font = 2) {
  graphics::rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2,
                 col = fill, border = border, lwd = 1.2, lty = lty)
  if (is.null(body)) {
    graphics::text(x, y, title, cex = title_cex, font = font, col = cols$ink)
  } else {
    graphics::text(x, y + h * 0.19, title, cex = title_cex, font = font,
                   col = cols$ink)
    graphics::text(x, y - h * 0.18, wrap(body), cex = body_cex,
                   col = cols$muted, lines = 0.92)
  }
}

arrow <- function(x0, y0, x1, y1, lty = 1, col = cols$ink) {
  graphics::arrows(x0, y0, x1, y1, length = 0.08, angle = 20,
                   lwd = 1.2, col = col, lty = lty)
}

elbow_arrow <- function(x, y, lty = 1, col = cols$ink) {
  n <- length(x)
  if (n < 2L || length(y) != n) stop("x and y must have the same length >= 2")
  if (n > 2L) {
    for (i in seq_len(n - 2L)) {
      graphics::segments(x[i], y[i], x[i + 1L], y[i + 1L],
                         lwd = 1.2, col = col, lty = lty)
    }
  }
  graphics::arrows(x[n - 1L], y[n - 1L], x[n], y[n],
                   length = 0.08, angle = 20, lwd = 1.2, col = col, lty = lty)
}

stage <- function(x, label, fill, width = NULL, y = 6.945) {
  if (is.null(width)) width <- max(0.64, 0.13 * nchar(label))
  graphics::rect(x - width / 2, y - 0.135, x + width / 2, y + 0.135,
                 col = fill, border = NA)
  graphics::text(x, y, toupper(label), cex = 0.55, font = 2,
                 col = cols$ink)
}

draw_map <- function() {
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i",
                family = "sans")
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 14), ylim = c(0, 8), asp = NA)

  graphics::rect(0, 0, 14, 8, col = "white", border = NA)
  graphics::text(0.55, 7.58, "bayesTLS function map", adj = 0,
                 cex = 1.15, font = 2, col = cols$ink)
  graphics::text(0.55, 7.25,
                 "Core user workflow only; report/access and prediction are parallel downstream uses.",
                 adj = 0, cex = 0.63, col = cols$muted)

  stage(1.25, "data", cols$data)
  stage(3.25, "fit", cols$fit)
  stage(5.65, "derive", cols$extract)
  stage(8.1, "report / access", cols$output, width = 1.58)
  stage(8.1, "predict", cols$predict, width = 0.88, y = 4.12)
  stage(11.45, "plot", cols$plot, width = 0.70)

  box(0.45, 5.2, 0.8, 0.72, "raw data", "counts or proportions",
      fill = "white", border = cols$grid, title_cex = 0.7, body_cex = 0.52,
      font = 1)
  box(1.85, 5.2, 1.65, 0.92, "standardize_data()",
      "schema + metadata", fill = cols$data)

  box(3.55, 5.2, 1.65, 0.92, "fit_4pl()",
      "joint Bayesian 4PL workflow", fill = cols$fit)
  box(3.55, 3.98, 1.65, 0.56, "get_4pl_est()",
      "4PL parameter summaries/draws", fill = cols$output,
      title_cex = 0.58, body_cex = 0.41)
  box(3.55, 3.20, 1.65, 0.58, "diagnostics",
      "diagnose_tdt_fit()\nbayes_R2_tls()", fill = cols$optional,
      border = cols$grid, title_cex = 0.55, body_cex = 0.40, lty = 2)
  box(3.55, 6.18, 1.65, 0.62, "model spec helpers",
      "make_4pl_formula()\nmake_4pl_priors()", fill = cols$optional,
      border = cols$grid, title_cex = 0.57, body_cex = 0.45, lty = 2)

  box(5.6, 5.75, 1.75, 0.86, "tls()",
      "general TLS extractor\nbayes_tls or brmsfit", fill = cols$extract)
  box(5.6, 4.78, 1.75, 0.86, "extract_tdt()",
      "workflow bundle\nz, CTmax, Tcrit, LT curve", fill = cols$extract)
  box(5.6, 3.55, 2.05, 0.74, "advanced primitives",
      "derive_z()\nderive_tdt_curve()\nderive_temperature_for_duration()",
      fill = cols$optional, border = cols$grid, title_cex = 0.58,
      body_cex = 0.45, lty = 2)
  box(5.6, 6.65, 1.45, 0.54, "manual brmsfit",
      "optional input to tls()", fill = "white", border = cols$grid,
      title_cex = 0.55, body_cex = 0.43, lty = 2, font = 1)

  box(8.1, 5.85, 2.2, 0.74, "get_tls_est()",
      "tls object -> TLS summaries/draws", fill = cols$output,
      title_cex = 0.64, body_cex = 0.47)
  box(8.1, 4.85, 2.2, 0.74, "extract_tdt accessors",
      "get_tls_summary()\nget_tls_draws()", fill = cols$output,
      title_cex = 0.55, body_cex = 0.45)

  box(8.1, 3.55, 2.05, 0.76, "predict_survival_curves()",
      "static temp x duration grid", fill = cols$predict)
  box(8.1, 2.65, 2.05, 0.76, "predict_heat_injury()",
      "dynamic traces + optional repair", fill = cols$predict)
  box(8.1, 1.90, 2.05, 0.50, "trace and repair helpers",
      "make_temperature_scenarios()\nrepair_rate_schoolfield()",
      fill = cols$optional, border = cols$grid, title_cex = 0.47,
      body_cex = 0.37, lty = 2)

  box(11.45, 4.55, 2.15, 1.18, "plot_*() family",
      "visualise TLS summaries,\nsurvival grids, landscapes,\nHI trajectories",
      fill = cols$plot)

  box(3.0, 0.86, 1.55, 0.52, "ts_stage1()",
      "classical LT50 by temp", fill = cols$optional, border = cols$grid,
      title_cex = 0.55, body_cex = 0.41, lty = 2)
  box(5.0, 0.86, 1.55, 0.52, "ts_stage2()",
      "OLS TDT line", fill = cols$optional, border = cols$grid,
      title_cex = 0.55, body_cex = 0.41, lty = 2)
  box(7.0, 0.86, 1.55, 0.52, "ts_ci() / ts_curve()",
      "uncertainty + line", fill = cols$optional, border = cols$grid,
      title_cex = 0.50, body_cex = 0.39, lty = 2)
  graphics::text(0.55, 0.86, "Comparison path\nnot the core Bayesian workflow",
                 adj = 0, cex = 0.56, col = cols$muted, lines = 0.95)

  arrow(0.85, 5.2, 1.02, 5.2, col = cols$grid)
  arrow(2.68, 5.2, 2.72, 5.2)
  arrow(3.55, 5.87, 3.55, 5.68, lty = 2, col = cols$muted)
  arrow(3.55, 4.74, 3.55, 4.27)
  arrow(4.38, 5.2, 4.70, 5.55)
  arrow(4.38, 5.2, 4.70, 4.86)
  arrow(5.6, 6.38, 5.6, 6.08, lty = 2, col = cols$muted)
  arrow(5.6, 4.35, 5.6, 3.93, lty = 2, col = cols$muted)
  arrow(6.48, 5.72, 7.00, 5.85)
  arrow(6.48, 4.78, 7.00, 4.85)
  elbow_arrow(c(4.38, 4.55, 7.00, 7.08), c(4.78, 2.98, 2.98, 3.55))
  elbow_arrow(c(4.38, 4.65, 7.00, 7.08), c(4.62, 2.22, 2.22, 2.65))
  arrow(9.20, 5.85, 10.38, 4.78)
  arrow(9.20, 4.85, 10.38, 4.55)
  arrow(9.13, 3.55, 10.38, 4.30)
  arrow(9.13, 2.65, 10.38, 4.05)
  arrow(8.1, 2.15, 8.1, 2.27, lty = 2, col = cols$muted)
  arrow(3.78, 0.86, 4.22, 0.86, lty = 2, col = cols$muted)
  arrow(5.78, 0.86, 6.22, 0.86, lty = 2, col = cols$muted)

  graphics::text(13.45, 0.45,
                 "Solid arrows = main workflow   Dashed = optional/advanced",
                 adj = 1, cex = 0.52, col = cols$muted)
}

svg_file <- file.path(out_dir, "bayesTLS_function_map.svg")
png_file <- file.path(out_dir, "bayesTLS_function_map.png")

grDevices::svg(svg_file, width = 14, height = 8, pointsize = 12)
draw_map()
grDevices::dev.off()

png_args <- list(filename = png_file, width = 2800, height = 1600, res = 200)
if (capabilities("cairo")) png_args$type <- "cairo"
do.call(grDevices::png, png_args)
draw_map()
grDevices::dev.off()

message("Wrote: ", svg_file)
message("Wrote: ", png_file)
