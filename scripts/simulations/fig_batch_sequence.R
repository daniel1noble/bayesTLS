# Experimental sequence behind day (batch) effects: the same 150 cups (5 temperatures x 6 exposures x 5 cups)
# run as A) one temperature per day, B) every temperature every day, or C) each temperature on two days.
# A day's shift (fill) is shared by every cup it holds. Companion to fig_batch_design.R; the simulation
# is sim_batch.R.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
temps   <- c(30, 32, 34, 36, 38)
shift_T <- c(0.06, -0.03, -0.10, 0.08, 0.02)      # illustrative day shifts (log10 min), days 1-5
shift_2 <- c(-0.02, 0.05, 0.03, -0.06, 0.07)      # days 6-10 in C
labA <- "A) One temperature per day (each day runs one temperature, 5 cups at each of 6 exposures)"
labB <- "B) Every temperature every day (1 cup per exposure each day, so 5 days give 5 cups)"
labC <- "C) Each temperature on two days (one temperature a day, its 5 cups per exposure split 3 + 2)"
labs <- c(labA, labB, labC)
runs <- bind_rows(
  tibble(design = labA, day = 1:5, temp = temps, shift = shift_T, cups = 5),
  expand_grid(design = labB, day = 1:5, temp = temps) |> mutate(shift = shift_T[day], cups = 1),
  tibble(design = labC, day = 1:10, temp = rep(temps, 2), shift = c(shift_T, shift_2), cups = rep(c(3, 2), each = 5)))
cups <- runs |>                                     # dots = cups: columns = 6 exposures, rows = replicate cups
  uncount(cups, .id = "cup", .remove = FALSE) |>
  expand_grid(dur = 1:6) |>
  mutate(dy = case_when(design == labA ~ (cup - 3) * 0.105 - 0.13,
                        design == labB ~ -0.08,
                        TRUE ~ (cup - (cups + 1) / 2) * 0.16),
         x = day + (dur - 3.5) * 0.12, y = match(temp, temps) + dy)
durs <- c(1, 5, 15, 45, 135, 405)                  # exposure times (min), labelled on one tile in A and B
keys <- bind_rows(tibble(design = labA, x = 1 + (1:6 - 3.5) * 0.12, y = 1.23, lab = as.character(durs)),
                  tibble(design = labA, x = 1, y = 1.37, lab = "exposure time (min)"),
                  tibble(design = labB, x = 1 + (1:6 - 3.5) * 0.12, y = 5.10, lab = as.character(durs)),
                  tibble(design = labB, x = 1, y = 5.30, lab = "exposure time (min)"))
runs <- mutate(runs, design = factor(design, labs))
cups <- mutate(cups, design = factor(design, labs))
keys <- mutate(keys, design = factor(design, labs))
p <- ggplot() +
  geom_tile(data = runs, aes(day, match(temp, temps), fill = shift), width = 0.9, height = 0.9,
            colour = "grey40", linewidth = 0.3) +
  geom_point(data = cups, aes(x, y), shape = 21, fill = "white", colour = "grey15", size = 1.1, stroke = 0.25) +
  geom_text(data = keys, aes(x, y, label = lab), size = 2, colour = "grey20") +
  facet_wrap(~ design, ncol = 1, scales = "free_x") +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-0.12, 0.12),
                       name = "Day's shift\n(log10 min)", breaks = c(-0.1, 0, 0.1),
                       labels = c("-0.1 died sooner", "0", "+0.1 lasted longer")) +
  scale_x_continuous(breaks = 1:10, labels = paste("Day", 1:10), expand = c(0.02, 0.02)) +
  scale_y_continuous(breaks = 1:5, labels = paste0(temps, " °C"), expand = c(0.02, 0.02)) +
  labs(x = NULL, y = "Bath temperature",
       caption = paste("All three designs test the same 150 cups. Each tile is one run (one bath on one day) and each white dot is one cup.",
                       "Across a tile are the 6 exposure times (1 to 405 min, labelled on the Day 1 tile in A and B), and down are replicate",
                       "cups at the same exposure.",
                       "A. Each temperature has its own day, so day and temperature are confounded. A day's shift hits every cup at that",
                       "temperature, so temperatures scatter off the TDT line (design figure b). The model needs (1 | day), which here is the",
                       "same grouping as temperature, and the day SD rests on how 5 temperatures scatter around a straight line.",
                       "B. Every temperature runs every day, so a day's shift moves the whole line up or down. The model needs (1 | day).",
                       "Because each day holds every temperature, day-wide effects (cohort, handler) can be told apart from temperature.",
                       "C. Each temperature runs on two separate days (3 cups per exposure on the first, 2 on the second). The model needs",
                       "(1 | day), which can now compare the two days at each temperature, so it measures day-to-day variation directly",
                       "rather than from the line alone. Each temperature's average also mixes two days' shifts.", sep = "\n")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold", hjust = 0),
        plot.caption = element_text(hjust = 0, colour = "grey30"), legend.position = "right")
ggsave(here::here("output", "sim_batch", "sequence_figure.png"), p, width = 9.5, height = 11.5, dpi = 150, bg = "white")
cat("saved\n")
