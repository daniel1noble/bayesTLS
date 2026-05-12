# Plotting helpers for the TDT function library. One plot function per
# prediction/extraction output, plus a shared theme. Each function returns a
# ggplot object the caller can extend with `+ <geom>` if needed.

#' Project ggplot theme
#'
#' A minimal classic theme with bold axis titles, used by all plotting
#' helpers in this library. Caller can override with `+ theme(...)`.
#'
#' @param base_size Base font size, passed to `theme_classic()`. Default 13.
#' @return A ggplot2 theme object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, hp)) + geom_point() + theme_tdt()
#' @export
theme_tdt <- function(base_size = 13) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.title       = ggplot2::element_text(face = "bold"),
      panel.border     = ggplot2::element_rect(colour = "black",
                                               fill = NA, linewidth = 0.6),
      legend.position  = "right"
    )
}

#' Plot posterior survival curves
#'
#' Plots posterior median ± 95% CrI of survival as a function of duration, one
#' coloured line per assay temperature. Optionally overlays observed survival
#' proportions per replicate.
#'
#' @param pred     Output of [predict_survival_curves()].
#' @param observed Optional standardised data tibble (from [standardize_data()])
#'                 to overlay as points (jittered, transparent).
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' pred <- predict_survival_curves(wf, temps = c(30, 32, 34, 36))
#' plot_survival_curves(pred, observed = wf$data)
#' }
#' @export
plot_survival_curves <- function(pred, observed = NULL) {
  df <- pred$summary
  p <- ggplot2::ggplot(df, ggplot2::aes(x = duration, y = survival_median,
                                        colour = factor(temp),
                                        fill   = factor(temp))) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = survival_lower,
                                      ymax = survival_upper),
                         alpha = 0.18, colour = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Exposure duration", y = "Survival probability",
                  colour = "Temperature (°C)", fill = "Temperature (°C)") +
    theme_tdt()

  if (!is.null(observed)) {
    p <- p + ggplot2::geom_jitter(
      data = observed,
      ggplot2::aes(x = duration, y = survival, colour = factor(temp)),
      width = 0.02, height = 0, alpha = 0.35, size = 0.8,
      inherit.aes = FALSE
    )
  }
  p
}

#' Plot a posterior LT_x curve
#'
#' Plots posterior median ± 95% CrI of the duration to reach `target_surv`
#' across temperature, with a log-scaled y axis (the classical TDT view).
#'
#' @param ltx Output of [derive_ltx_curve()].
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' ltx <- derive_ltx_curve(wf, temp_grid = seq(29, 37, by = 0.5))
#' plot_ltx_curve(ltx)
#' }
#' @export
plot_ltx_curve <- function(ltx) {
  target_lab <- paste0(round(100 * unique(ltx$summary$target_surv)),
                       "% survival")
  ggplot2::ggplot(ltx$summary,
                  ggplot2::aes(x = temp, y = duration_median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = duration_lower,
                                      ymax = duration_upper),
                         fill = "#146C7C", alpha = 0.2) +
    ggplot2::geom_line(colour = "#146C7C", linewidth = 1) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Assay temperature (°C)",
                  y = paste0("Time to ", target_lab, " (",
                             ltx$output_time_unit, ")")) +
    theme_tdt()
}

#' Plot the TDT survival landscape as a heatmap
#'
#' @param landscape Output of [derive_tdt_landscape()].
#' @param observed  Optional standardised data tibble to overlay.
#' @param contours  Numeric vector of survival levels to draw contours at.
#'                  Default `c(0.25, 0.5, 0.75)`.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' lsp <- derive_tdt_landscape(wf)
#' plot_tdt_landscape(lsp, observed = wf$data)
#' }
#' @export
plot_tdt_landscape <- function(landscape, observed = NULL,
                               contours = c(0.25, 0.5, 0.75)) {
  df <- landscape$summary
  p <- ggplot2::ggplot(df, ggplot2::aes(x = temp, y = duration,
                                        fill = survival_median)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::geom_contour(ggplot2::aes(z = survival_median),
                          breaks = contours,
                          colour = "white", alpha = 0.7) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), option = "inferno") +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Temperature (°C)",
                  y = "Exposure duration",
                  fill = "Survival") +
    theme_tdt()

  if (!is.null(observed)) {
    p <- p + ggplot2::geom_point(
      data = observed,
      ggplot2::aes(x = temp, y = duration, fill = survival),
      inherit.aes = FALSE, shape = 21, colour = "white",
      size = 2, stroke = 0.2
    )
  }
  p
}

#' Plot a posterior temperature density (e.g. CTmax or T_crit)
#'
#' Density plot of a per-draw temperature posterior, with a horizontal segment
#' marking the 95% CrI and a point at the median. Useful for visualising
#' [derive_temperature_for_duration()] output, or the `$CTmax` / `$T_crit`
#' elements of [extract_tdt()].
#'
#' @param temp_post  Output of [derive_temperature_for_duration()] OR a list
#'                   with `$draws` and `$summary` from [extract_tdt()].
#' @param truth      Optional numeric scalar: a true value to mark with a
#'                   dashed vertical line.
#' @param x_label    X-axis label. Default `"Temperature (°C)"`.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' et <- extract_tdt(wf)
#' plot_temperature_density(et$CTmax)
#' plot_temperature_density(et$T_crit)
#' }
#' @export
plot_temperature_density <- function(temp_post, truth = NULL,
                                     x_label = "Temperature (°C)") {
  draws <- temp_post$draws |> dplyr::filter(is.finite(temp))
  d     <- stats::density(draws$temp, na.rm = TRUE)
  y_ci  <- -max(d$y, na.rm = TRUE) * 0.1

  p <- ggplot2::ggplot(draws, ggplot2::aes(x = temp)) +
    ggplot2::geom_density(fill = "#146C7C", alpha = 0.35,
                          colour = "#0E4A55", linewidth = 0.4) +
    ggplot2::geom_segment(
      data = temp_post$summary,
      ggplot2::aes(x = temp_lower, xend = temp_upper,
                   y = y_ci, yend = y_ci),
      linewidth = 1, colour = "#2F2F2F", inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = temp_post$summary,
      ggplot2::aes(x = temp_median, y = y_ci),
      size = 2.5, colour = "#2F2F2F", inherit.aes = FALSE
    ) +
    ggplot2::labs(x = x_label, y = "Posterior density") +
    theme_tdt()

  if (!is.null(truth)) {
    p <- p + ggplot2::geom_vline(xintercept = truth,
                                  linetype = "dashed", colour = "grey40")
  }
  p
}

#' Plot the three reference temperature scenarios
#'
#' One panel per scenario stacked vertically, with optional horizontal dashed
#' line at the damage threshold `T_c`.
#'
#' @param scens Named list from [make_temperature_scenarios()].
#' @param T_c   Optional damage threshold to mark.
#' @return A ggplot object.
#' @examples
#' scens <- make_temperature_scenarios()
#' plot_temperature_scenarios(scens, T_c = 24)
#' @export
plot_temperature_scenarios <- function(scens, T_c = NULL) {
  df <- dplyr::bind_rows(
    lapply(names(scens), function(nm) {
      x <- scens[[nm]]; x$scenario <- nm; x
    })
  )
  df$scenario <- factor(df$scenario,
                        levels = c("flat", "single_spike", "multi_spike"))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = time_h, y = temp)) +
    ggplot2::geom_line(linewidth = 0.7, colour = "#146C7C") +
    ggplot2::facet_wrap(~ scenario, ncol = 1) +
    ggplot2::labs(x = "Time (hours)", y = "Temperature (°C)") +
    theme_tdt()

  if (!is.null(T_c)) {
    p <- p + ggplot2::geom_hline(yintercept = T_c,
                                  linetype = "dashed",
                                  colour = "grey50")
  }
  p
}

#' Plot HI and survival trajectories from `predict_heat_injury()`
#'
#' Two stacked panels: cumulative heat injury (%) on top, predicted survival
#' fraction on the bottom, both with 95% credible bands.
#'
#' @param hi             Output of [predict_heat_injury()].
#' @param lt50_threshold Numeric. Horizontal dashed line on the HI panel
#'                       marking the LT50 dose threshold. Default 100.
#' @return A ggplot object (combined via `patchwork`).
#' @examples
#' \dontrun{
#' hi <- predict_heat_injury(scens$single_spike, wf)
#' plot_heat_injury(hi)
#' }
#' @export
plot_heat_injury <- function(hi, lt50_threshold = 100) {
  df <- hi$summary

  p_hi <- ggplot2::ggplot(df, ggplot2::aes(x = time_h, y = hi_median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = hi_lower, ymax = hi_upper),
                         fill = "#C76D37", alpha = 0.2) +
    ggplot2::geom_line(colour = "#C76D37", linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = lt50_threshold,
                         linetype = "dashed", colour = "grey40") +
    ggplot2::labs(x = NULL, y = "Cumulative HI (%)") +
    theme_tdt()

  p_surv <- ggplot2::ggplot(df, ggplot2::aes(x = time_h, y = surv_median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = surv_lower, ymax = surv_upper),
                         fill = "#146C7C", alpha = 0.2) +
    ggplot2::geom_line(colour = "#146C7C", linewidth = 0.9) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Time (hours)", y = "Predicted survival") +
    theme_tdt()

  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p_hi, p_surv, ncol = 1)
  } else {
    list(hi = p_hi, survival = p_surv)
  }
}

#' Plot a Sharpe-Schoolfield repair TPC
#'
#' @param temp_grid   Numeric vector of temperatures (°C) to evaluate.
#' @param repair_pars Named list of Sharpe-Schoolfield parameters
#'                    (`TA, TAL, TAH, TL, TH, TREF, r_ref`).
#' @return A ggplot object.
#' @examples
#' rp <- list(TA = 14065, TAL = 50000, TAH = 120000,
#'            TL = 10.5 + 273.15, TH = 22.5 + 273.15,
#'            TREF = 17 + 273.15, r_ref = 0.01)
#' plot_repair_rate(seq(5, 30, length.out = 200), rp)
#' @export
plot_repair_rate <- function(temp_grid, repair_pars) {
  rates <- repair_rate_schoolfield(
    temp_celsius = temp_grid,
    TA = repair_pars$TA, TAL = repair_pars$TAL,
    TAH = repair_pars$TAH, TL = repair_pars$TL,
    TH = repair_pars$TH, TREF = repair_pars$TREF,
    r_ref = repair_pars$r_ref
  )
  df <- tibble::tibble(temp = temp_grid, repair_rate = rates)

  ggplot2::ggplot(df, ggplot2::aes(x = temp, y = repair_rate)) +
    ggplot2::geom_line(linewidth = 1, colour = "#4DAC26") +
    ggplot2::labs(x = "Temperature (°C)",
                  y = "Repair rate (per hour)") +
    theme_tdt()
}
