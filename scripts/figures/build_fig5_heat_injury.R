#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(bayesTLS)
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(ggtext)
  library(here)
  library(magick)
  library(patchwork)
  library(readr)
  library(scales)
})

out_fig <- here::here("output", "figs")
out_dat <- here::here("output", "data")
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dat, showWarnings = FALSE, recursive = TRUE)

objs_path <- here::here("output", "figs", "_fig_objs.rds")
if (!file.exists(objs_path)) {
  stop("Missing output/figs/_fig_objs.rds. Run the full manuscript figure build first.",
       call. = FALSE)
}
list2env(readRDS(objs_path), environment())

P <- function(x) here::here("pics", x)
svg_png <- local({
  cache <- new.env(parent = emptyenv())
  function(svg, width = 1600) {
    key <- paste0(svg, "@", width)
    if (is.null(cache[[key]])) {
      out <- tempfile(fileext = ".png")
      magick::image_write(magick::image_read_svg(P(svg), width = width), out)
      cache[[key]] <- out
    }
    cache[[key]]
  }
})

message(">> building Fig. 5 heat-injury data with r_ref = 0.005")

ND <- 150
set.seed(123)

orsted_micro <- utils::read.csv(here::here(
  "inst", "extdata", "orsted_2024",
  "orsted2024_nichemapr_rennes_2018_hourly.csv.gz"
))
fly_trace <- orsted_micro |>
  mutate(datetime_utc = as.POSIXct(datetime_utc, tz = "UTC")) |>
  filter(substr(datetime_utc, 1, 4) == "2018",
         substr(datetime_utc, 6, 7) %in% c("06", "07", "08"),
         is.finite(micro_temp_c)) |>
  arrange(datetime_utc) |>
  transmute(time = row_number() - 1, temp = micro_temp_c,
            datetime = datetime_utc)

Tc_F <- get_tls_est(tls(wf_dros_mort_F, target_surv = "relative", t_ref = 60,
                        lethal = TRUE, ndraws = 150), "summary", "Tcrit")$median
Tc_M <- get_tls_est(tls(wf_dros_mort_M, target_surv = "relative", t_ref = 60,
                        lethal = TRUE, ndraws = 150), "summary", "Tcrit")$median

fly_repair <- list(TA = 14065, TAL = 50000, TAH = 100000,
                   TL = 15 + 273.15, TH = 28 + 273.15,
                   TREF = 25 + 273.15, r_ref = 0.005)

agg_draws <- function(hi, ...) {
  hi$draws |>
    group_by(time) |>
    summarise(hi_med = median(hi), hi_se = sd(hi),
              surv_med = median(survival), surv_se = sd(survival),
              .groups = "drop") |>
    mutate(...)
}

fly_rows <- list()
for (sx in c("Female", "Male")) {
  wf <- if (sx == "Female") wf_dros_mort_F else wf_dros_mort_M
  Tc <- if (sx == "Female") Tc_F else Tc_M
  for (rep in c(FALSE, TRUE)) {
    hi <- predict_heat_injury(
      fly_trace[, c("time", "temp")], wf,
      target_surv = "relative", T_c = Tc, trace_unit = "hours",
      ndraws = ND, repair = rep,
      repair_pars = if (rep) fly_repair else NULL, save_draws = TRUE
    )
    fly_rows[[length(fly_rows) + 1]] <-
      agg_draws(hi, sex = sx, repair = ifelse(rep, "With repair", "No repair"))
  }
}
fdat <- bind_rows(fly_rows) |> left_join(fly_trace, by = "time")
fTc <- mean(c(Tc_F, Tc_M))

aphid_trace <- read_csv(
  here::here("inst", "extdata", "data_temp_trace_aphid_summer2016.csv"),
  show_col_types = FALSE
) |>
  filter(city == "Wuhan", as.integer(format(datetime, "%m")) <= 6) |>
  transmute(time = time_h, temp = temp_c,
            datetime = as.POSIXct(sub("T", " ", as.character(datetime)), tz = "UTC"))

aphid_repair <- list(TA = 14065, TAL = 50000, TAH = 100000,
                     TL = 15 + 273.15, TH = 35 + 273.15,
                     TREF = 25 + 273.15, r_ref = 0.005)
aphid_Tc <- aphid_tls$summary |>
  filter(quantity == "Tcrit") |>
  transmute(species = as.character(species), Tc = median)

aphid_rows <- list()
for (rep in c(FALSE, TRUE)) {
  hi <- predict_heat_injury(
    aphid_trace[, c("time", "temp")], wf_aphid,
    target_surv = "relative", trace_unit = "hours", by = "species",
    ndraws = ND, repair = rep,
    repair_pars = if (rep) aphid_repair else NULL, save_draws = TRUE
  )
  aphid_rows[[length(aphid_rows) + 1]] <- hi$draws |>
    group_by(species, time) |>
    summarise(hi_med = median(hi), hi_se = sd(hi),
              surv_med = median(survival), surv_se = sd(survival),
              .groups = "drop") |>
    mutate(repair = ifelse(rep, "With repair", "No repair"))
}
adat <- bind_rows(aphid_rows) |>
  mutate(species = as.character(species)) |>
  left_join(aphid_trace, by = "time")

saveRDS(
  list(fly = fdat, fly_trace = fly_trace, fly_Tc = fTc,
       aphid = adat, aphid_trace = aphid_trace, aphid_Tc = aphid_Tc,
       repair_pars = list(fly = fly_repair, aphid = aphid_repair)),
  here::here("output", "data", "fig5_hi_data.rds")
)

fdat$repair <- factor(fdat$repair, levels = c("No repair", "With repair"))
adat$repair <- factor(adat$repair, levels = c("No repair", "With repair"))

sp_pal <- c(M_dirhodum = "#E69F00", S_avenae = "#D55E00", R_padi = "#CC79A7")
sp_lab <- c(M_dirhodum = "M. dirhodum", S_avenae = "S. avenae", R_padi = "R. padi")
sex_pal <- c(Female = "#009E73", Male = "#7B3294")
tcrit_col <- "#7B3294"

f6_thm <- theme_classic(base_size = 13) +
  theme(panel.grid.major.y = element_line(colour = "grey93"),
        plot.margin = margin(2, 7, 2, 3),
        legend.position = "none",
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))
f6_xsc <- function() scale_x_datetime(date_breaks = "1 month", date_labels = "%b")
f6_tg <- function(p, lab) {
  p + labs(tag = lab) +
    theme(plot.tag = element_text(face = "bold", size = 12),
          plot.tag.position = c(0.015, 0.98))
}

f6_temp1 <- function(tr, Tc, title, ylab) {
  ggplot(tr, aes(datetime, temp)) +
    geom_line(colour = "grey45", linewidth = 0.3) +
    geom_hline(yintercept = Tc, linetype = "dotted",
               colour = tcrit_col, linewidth = 0.9) +
    annotate("text", x = min(tr$datetime), y = Tc, label = "T[crit]",
             hjust = -0.1, vjust = -0.4, size = 3.7, colour = tcrit_col,
             fontface = "bold", parse = TRUE) +
    labs(y = ylab, x = NULL, title = title) +
    f6_xsc() + f6_thm +
    theme(axis.text.x = element_blank(),
          plot.title = element_markdown(face = "bold", size = 13, hjust = 0.5))
}

f6_tempA <- function(tr, tc_df, title, ylab) {
  ggplot(tr, aes(datetime, temp)) +
    geom_line(colour = "grey45", linewidth = 0.3) +
    geom_hline(data = tc_df, aes(yintercept = Tc, colour = species),
               linetype = "dotted", linewidth = 0.9, show.legend = FALSE) +
    annotate("text", x = min(tr$datetime), y = max(tc_df$Tc),
             label = "T[crit]", hjust = -0.1, vjust = -0.5,
             size = 3.7, colour = "grey15", fontface = "bold", parse = TRUE) +
    scale_colour_manual(values = sp_pal) +
    labs(y = ylab, x = NULL, title = title) +
    f6_xsc() + f6_thm +
    theme(axis.text.x = element_blank(),
          plot.title = element_markdown(face = "bold", size = 13, hjust = 0.5))
}

f6_hiF <- function(dd) {
  ggplot(dd, aes(datetime, hi_med, colour = sex, fill = sex)) +
    geom_ribbon(aes(ymin = pmax(0, hi_med - hi_se), ymax = hi_med + hi_se,
                    group = interaction(sex, repair)),
                colour = NA, alpha = 0.13) +
    geom_hline(yintercept = 100, linetype = "dotted", colour = "grey55") +
    geom_line(aes(linetype = repair), linewidth = 0.55) +
    scale_colour_manual(values = sex_pal) +
    scale_fill_manual(values = sex_pal, guide = "none") +
    scale_linetype_manual(values = c("No repair" = "solid", "With repair" = "22")) +
    labs(y = "Heat injury (%)", x = NULL) +
    f6_xsc() + f6_thm + theme(axis.text.x = element_blank())
}

f6_survF <- function(dd) {
  ggplot(dd, aes(datetime, surv_med, colour = sex, fill = sex)) +
    geom_ribbon(aes(ymin = pmax(0, surv_med - surv_se),
                    ymax = pmin(1, surv_med + surv_se),
                    group = interaction(sex, repair)),
                colour = NA, alpha = 0.13) +
    geom_line(aes(linetype = repair), linewidth = 0.55) +
    scale_colour_manual(values = sex_pal) +
    scale_fill_manual(values = sex_pal, guide = "none") +
    scale_linetype_manual(values = c("No repair" = "solid", "With repair" = "22")) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    labs(y = "Survival", x = NULL) +
    f6_xsc() + f6_thm
}

f6_hiA <- function(dd) {
  ggplot(dd, aes(datetime, hi_med, colour = species, fill = species)) +
    geom_ribbon(aes(ymin = pmax(0, hi_med - hi_se), ymax = hi_med + hi_se,
                    group = interaction(species, repair)),
                colour = NA, alpha = 0.12) +
    geom_hline(yintercept = 100, linetype = "dotted", colour = "grey55") +
    geom_line(aes(linetype = repair), linewidth = 0.55) +
    scale_colour_manual(values = sp_pal, labels = sp_lab) +
    scale_fill_manual(values = sp_pal, guide = "none") +
    scale_linetype_manual(values = c("No repair" = "solid", "With repair" = "22")) +
    labs(y = "Heat injury (%)", x = NULL) +
    f6_xsc() + f6_thm + theme(axis.text.x = element_blank())
}

f6_survA <- function(dd) {
  ggplot(dd, aes(datetime, surv_med, colour = species, fill = species)) +
    geom_ribbon(aes(ymin = pmax(0, surv_med - surv_se),
                    ymax = pmin(1, surv_med + surv_se),
                    group = interaction(species, repair)),
                colour = NA, alpha = 0.12) +
    geom_line(aes(linetype = repair), linewidth = 0.55) +
    scale_colour_manual(values = sp_pal, labels = sp_lab) +
    scale_fill_manual(values = sp_pal, guide = "none") +
    scale_linetype_manual(values = c("No repair" = "solid", "With repair" = "22")) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    labs(y = "Survival", x = NULL) +
    f6_xsc() + f6_thm
}

f6_ins <- function(p, f, l, b, r, t) {
  p + inset_element(ggdraw() + draw_image(svg_png(f)),
                    left = l, bottom = b, right = r, top = t, align_to = "panel")
}
f6_ins_png <- function(p, png, l, b, r, t) {
  p + inset_element(ggdraw() + draw_image(P(png)),
                    left = l, bottom = b, right = r, top = t, align_to = "panel")
}

fly_hi <- f6_ins(f6_tg(f6_hiF(fdat), "b)"), "dros_male_v2.svg",
                 .02, .55, .34, .98)
aphid_hi <- f6_ins_png(f6_tg(f6_hiA(adat), "e)"), "aphid_Mdirhodum_V2.png",
                       .02, .46, .44, .99)
col_f <- f6_tg(f6_temp1(fly_trace, fTc, "Vinegar fly (*D. suzukii*)",
                        "Temperature (deg C)"), "a)") /
  fly_hi / f6_tg(f6_survF(fdat), "c)")
col_a <- f6_tg(f6_tempA(aphid_trace, aphid_Tc, "Cereal aphids (Wuhan)",
                        "Temperature (deg C)"), "d)") /
  aphid_hi / f6_tg(f6_survA(adat), "f)")

mk_leg <- function(p, which) {
  g <- p + theme(legend.position = "bottom", legend.title = element_blank())
  if (which == "colour") {
    g <- g + guides(fill = "none", linetype = "none", colour = guide_legend(order = 1))
  } else {
    g <- g + guides(fill = "none", colour = "none", linetype = guide_legend(order = 1))
  }
  cowplot::get_legend(g)
}

fig5_final <- cowplot::plot_grid(
  (col_f | col_a) + plot_layout(heights = c(0.5, 1, 1)),
  cowplot::plot_grid(mk_leg(f6_hiF(fdat), "colour"),
                     mk_leg(f6_hiA(adat), "colour"),
                     mk_leg(f6_hiA(adat), "linetype"),
                     ncol = 3, rel_widths = c(1, 1.3, 1)),
  ncol = 1, rel_heights = c(1, 0.09)
)

ggsave(here::here("output", "figs", "fig5_heat_injury.png"), fig5_final,
       width = 9.8, height = 8.0, dpi = 300, bg = "white")

message(">> Fig. 5 heat-injury figure rebuilt")
message(">> saved output/figs/fig5_heat_injury.png")
message(">> saved output/data/fig5_hi_data.rds")
