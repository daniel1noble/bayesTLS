# Design figure for the non-independence simulation (sim_batch.R): the true time to 50% survival
# against temperature. Each day shifts every cup it holds by the same amount. The shifts drawn are the
# illustrative ones of the sequence figure (fig_batch_sequence.R), one per day.
source(here::here("scripts", "simulations", "sim_functions.R"))
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
tr    <- sim_truth("baseline", u_0 = 0.999, ell_0 = 0.001)              # the truth of sim_batch.R
temps <- design_spec("full")$temps
line  <- function(T) tr$m_beta0 + tr$m_beta1 * (T - tr$T_bar) +          # true log10 time to 50% survival
  log((tr$u - 0.5) / (0.5 - tr$ell)) / tr$k                               # (midpoint + offset to 50%)
shift <- c(0.06, -0.03, -0.10, 0.08, 0.02)                                 # one per day (log10 min)
shift2 <- c(-0.02, 0.05, 0.03, -0.06, 0.07)                                # second day per temperature in C
labs3 <- c("a) Truth: one straight TDT line", "b) A: one temperature per day", "c) B: every temperature every day",
           "d) C: each temperature on two days")
pts <- bind_rows(
  tibble(panel = labs3[1], T = temps, mid = line(temps)),
  tibble(panel = labs3[2], day = 1:5, T = temps, shift = shift) |> mutate(mid = line(T) + shift),
  tidyr::crossing(day = 1:5, T = temps) |> mutate(panel = labs3[3], shift = shift[day], mid = line(T) + shift),
  tibble(panel = labs3[4], day = 1:10, T = rep(temps, 2) + rep(c(-0.25, 0.25), each = 5),
         shift = c(shift, shift2)) |> mutate(mid = line(T) + shift)) |>
  mutate(panel = factor(panel, labs3))
ref  <- tidyr::crossing(panel = factor(labs3, labs3), T = c(29.5, 38.5)) |> mutate(mid = line(T))
days <- tidyr::crossing(day = 1:5, T = c(29.5, 38.5)) |>
  mutate(panel = factor(labs3[3], labs3), shift = shift[day], mid = line(T) + shift)
p <- ggplot(pts, aes(T, 10^mid)) +
  geom_line(data = ref, colour = "grey45", linewidth = 0.9) +
  geom_line(data = days, aes(group = day, colour = shift), linewidth = 0.6) +
  geom_segment(data = filter(pts, panel %in% labs3[c(2, 4)]),
               aes(xend = T, y = 10^line(T), yend = 10^mid, colour = shift), linewidth = 0.6) +
  geom_point(data = filter(pts, panel == labs3[1]), size = 2.4) +
  geom_point(data = filter(pts, panel != labs3[1]), aes(fill = shift), shape = 21, colour = "grey20", size = 2.4) +
  facet_wrap(~ panel, nrow = 2) +
  scale_colour_gradient2(low = "#B2182B", mid = "grey80", high = "#2166AC", midpoint = 0,
                         limits = c(-0.12, 0.12), guide = "none") +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-0.12, 0.12),
                       name = "Day's shift\n(log10 min)", breaks = c(-0.1, 0, 0.1),
                       labels = c("-0.1 died sooner", "0", "+0.1 lasted longer")) +
  scale_y_log10(breaks = c(5, 10, 20, 50, 100, 200, 500)) +
  scale_x_continuous(breaks = temps) +
  labs(x = "Assay temperature (°C)", y = "True time to 50% survival (min, log scale)",
       caption = paste(sprintf("Grey line = the straight TDT line every method aims at (z = %.1f °C, CTmax 1 h = %.1f °C).",
                               tr$z_true, tr$CTmax_1hr_true),
                       "Each day shifts every cup it holds by the same amount, drawn anew for each simulated dataset (SD 0.08 log10 min,",
                       "x1.2 in time). The shifts shown are those of the sequence figure.",
                       "b) A runs one temperature per day, so the shifts scatter the temperatures off the line and change the slope (z).",
                       "c) B runs every temperature every day, so a day's shift moves its whole line. The slope stays the same and only the",
                       "height of the line (CTmax) moves.",
                       "d) C runs each temperature on two days, so each temperature has two independent shifts. Their difference shows the",
                       "model how large day effects are, and their average sits closer to the line.", sep = "\n")) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0),
        plot.caption = element_text(hjust = 0, colour = "grey30"))
ggsave(here::here("output", "sim_batch", "design_figure.png"), p, width = 10, height = 8.6, dpi = 150, bg = "white")
cat("saved\n")
