
# =========================
# 1) Load + process raw data
# =========================
d <- read_excel("data/data_zebrafish_TDT.xlsx", sheet = "LETHAL_TDT")

d <- d %>%
  filter(Exclude == "no") %>%
  mutate(
    across(
      c(
        N_total,
        Temperature_assay,
        Duration_exposure_hours,
        N_malformed,
        matches("^(Hatching|Mortality)_day_\\d+_(morning|afternoon)$")
      ),
      ~ readr::parse_number(as.character(.))
    )
  ) %>%
  mutate(
    n_total = N_total,
    n_dead = rowSums(across(matches("^Mortality_day_\\d+_(morning|afternoon)$")), na.rm = TRUE),
    n_surv = n_total - n_dead,
    survival = n_surv / n_total,
    malformation_rate = N_malformed / pmax(n_surv, 1),
    assay_temp = round(Temperature_assay, 1),
    duration = as.numeric(Duration_exposure_hours),
    logd = log10(duration),
    life_stage = factor(Life_stage, levels = c("larvae", "old_embryos", "young_embryos")),
    date_experiment = factor(Date_experiment),
    sample = factor(Sample)
  ) %>%
  filter(duration > 0) # remove controls

d$obs <- factor(seq_len(nrow(d)))

temp_mean <- mean(d$assay_temp, na.rm = TRUE)

d <- d %>%
  mutate(
    temp_c = assay_temp - temp_mean # center
  )

# =========================
# 2) Numeric constants for priors / init
# =========================
mid_start <- median(d$logd, na.rm = TRUE)
lraw_mean <- qlogis(0.02 / 0.49)
uraw_mean <- qlogis((0.98 - 0.51) / 0.49)
logk_mean <- log(2)

num_to_string <- function(x) {
  formatC(x, digits = 16, format = "fg", flag = "#")
}

mid_start_s <- num_to_string(mid_start)
lraw_mean_s <- num_to_string(lraw_mean)
uraw_mean_s <- num_to_string(uraw_mean)
logk_mean_s <- num_to_string(logk_mean)

# =========================
# 3) 4-parameter logistic model over duration
# =========================
# Probability-scale curve:
#
# p = low + (up - low) / (1 + exp(exp(logk) * (logd - mid)))
#
# low = inv_logit(lraw) * 0.49
# up  = 0.51 + inv_logit(uraw) * 0.49
#
# mid varies with temperature and life stage:
# mid ~ temp_c * life_stage + (1 | date_experiment)
#
# z is an obs-level random effect

f_tdt_4pl <- bf(
  n_surv | trials(n_total) ~ eta,
  nlf(
    eta ~ logit(
      (inv_logit(lraw) * 0.49) +
        (
          (0.51 + inv_logit(uraw) * 0.49) -
            (inv_logit(lraw) * 0.49)
        ) / (1 + exp(exp(logk) * (logd - mid)))
    ) + z
  ),
  lf(lraw ~ 0 + life_stage),
  lf(uraw ~ 0 + life_stage),
  lf(logk ~ 0 + life_stage),
  lf(mid ~ temp_c * life_stage + (1 | date_experiment)),
  lf(z ~ 0 + (1 | obs)),
  nl = TRUE
)

# =========================
# 4) Priors
# =========================
priors_tdt_4pl <- c(
  set_prior(
    paste0("normal(", lraw_mean_s, ", 1)"),
    class = "b", nlpar = "lraw"
  ),
  set_prior(
    paste0("normal(", uraw_mean_s, ", 1)"),
    class = "b", nlpar = "uraw"
  ),
  set_prior(
    paste0("normal(", logk_mean_s, ", 1)"),
    class = "b", nlpar = "logk"
  ),
  set_prior(
    paste0("normal(", mid_start_s, ", 1.5)"),
    class = "b", nlpar = "mid", coef = "Intercept"
  ),
  set_prior(
    "normal(0, 0.6)",
    class = "b", nlpar = "mid"
  ),
  set_prior(
    "exponential(2)",
    class = "sd", nlpar = "mid", group = "date_experiment"
  ),
  set_prior(
    "exponential(2)",
    class = "sd", nlpar = "z", group = "obs"
  )
)

# =========================
# 5) Initial values
# =========================
init_tdt_4pl <- function() {
  list(
    b_lraw_life_stagelarvae = lraw_mean,
    b_lraw_life_stageold_embryos = lraw_mean,
    b_lraw_life_stageyoung_embryos = lraw_mean,
    b_uraw_life_stagelarvae = uraw_mean,
    b_uraw_life_stageold_embryos = uraw_mean,
    b_uraw_life_stageyoung_embryos = uraw_mean,
    b_logk_life_stagelarvae = logk_mean,
    b_logk_life_stageold_embryos = logk_mean,
    b_logk_life_stageyoung_embryos = logk_mean,
    b_mid_Intercept = mid_start
  )
}

# =========================
# 6) Fit model
# =========================
fit_bayesian_4pl <- brm(
  formula = f_tdt_4pl,
  data = d,
  family = binomial(link = "logit"),
  prior = priors_tdt_4pl,
  chains = 4,
  cores = 4,
  iter = 6000,
  warmup = 3000,
  init = init_tdt_4pl,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.995, max_treedepth = 15),
  seed = 123
)

print(fit_bayesian_4pl)

# =========================
# 7) Comparison with no obs overdispersion term
# =========================
f_tdt_4pl_noobs <- bf(
  n_surv | trials(n_total) ~ eta,
  nlf(
    eta ~ logit(
      (inv_logit(lraw) * 0.49) +
        (
          (0.51 + inv_logit(uraw) * 0.49) -
            (inv_logit(lraw) * 0.49)
        ) / (1 + exp(exp(logk) * (logd - mid)))
    )
  ),
  lf(lraw ~ 0 + life_stage),
  lf(uraw ~ 0 + life_stage),
  lf(logk ~ 0 + life_stage),
  lf(mid ~ temp_c * life_stage + (1 | date_experiment)),
  nl = TRUE
)

priors_tdt_4pl_noobs <- c(
  set_prior(
    paste0("normal(", lraw_mean_s, ", 1)"),
    class = "b", nlpar = "lraw"
  ),
  set_prior(
    paste0("normal(", uraw_mean_s, ", 1)"),
    class = "b", nlpar = "uraw"
  ),
  set_prior(
    paste0("normal(", logk_mean_s, ", 1)"),
    class = "b", nlpar = "logk"
  ),
  set_prior(
    paste0("normal(", mid_start_s, ", 1.5)"),
    class = "Intercept", nlpar = "mid"
  ),
  set_prior(
    "normal(0, 0.6)",
    class = "b", nlpar = "mid"
  ),
  set_prior(
    "exponential(2)",
    class = "sd", nlpar = "mid", group = "date_experiment"
  )
)

fit_bayesian_4pl_noobs <- brm(
  formula = f_tdt_4pl_noobs,
  data = d,
  family = binomial(link = "logit"),
  prior = priors_tdt_4pl_noobs,
  chains = 4,
  cores = 4,
  iter = 6000,
  warmup = 3000,
  init = init_tdt_4pl,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.995, max_treedepth = 15),
  seed = 123
)

loo_full <- loo(fit_bayesian_4pl)
loo_noobs <- loo(fit_bayesian_4pl_noobs)
loo_compare(loo_full, loo_noobs) # Keep the observation-level random effect for overdispersion

# =========================
# 8) Survival curves
# =========================
temps_plot <- seq(38, 42, by = 0.5)

dur_min <- min(d$duration, na.rm = TRUE)
dur_max <- max(d$duration, na.rm = TRUE)
dur_grid <- exp(seq(log(dur_min), log(dur_max), length.out = 250))

nd_surv <- expand.grid(
  assay_temp = temps_plot,
  duration = dur_grid,
  life_stage = levels(d$life_stage),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

nd_surv <- nd_surv %>%
  mutate(
    logd = log10(duration),
    temp_c = assay_temp - temp_mean,
    n_total = as.integer(round(median(d$n_total, na.rm = TRUE))),
    date_experiment = levels(d$date_experiment)[1],
    obs = levels(d$obs)[1],
    temp_f = factor(assay_temp, levels = sort(unique(assay_temp)))
  )

pr_surv <- posterior_linpred(
  fit_bayesian_4pl,
  newdata = nd_surv,
  re_formula = NA,
  transform = TRUE
)

nd_surv$pred_mean <- apply(pr_surv, 2, median, na.rm = TRUE)
nd_surv$pred_lo   <- apply(pr_surv, 2, quantile, probs = 0.025, na.rm = TRUE)
nd_surv$pred_hi   <- apply(pr_surv, 2, quantile, probs = 0.975, na.rm = TRUE)

p_surv_curves <- ggplot(
  nd_surv,
  aes(x = duration, y = pred_mean, colour = temp_f, fill = temp_f, group = temp_f)
) +
  geom_ribbon(
    aes(ymin = pred_lo, ymax = pred_hi),
    alpha = 0.15,
    colour = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 1.1) +
  facet_wrap(
    ~ life_stage,
    labeller = labeller(life_stage = life_stage_labels)
  ) +
  scale_colour_viridis_d(option = "plasma", direction = 1, begin = 0.05, end = 0.95) +
  scale_fill_viridis_d(option = "plasma", direction = 1, begin = 0.05, end = 0.95) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Exposure duration (hours)",
    y = "Predicted survival",
    colour = "Assay temperature (°C)"
  ) +
  base_theme +
  theme(
    legend.position = "right"
  )

print(p_surv_curves)


# =========================
# 4) Posterior draws helper functions
#    for 4PL TDT model
# =========================
dr <- posterior::as_draws_df(fit_bayesian_4pl)

get_par <- function(nm) {
  if (nm %in% names(dr)) {
    dr[[nm]]
  } else {
    rep(0, nrow(dr))
  }
}

get_mid_draws <- function(stage_name, temp_c_value) {
  if (stage_name == "larvae") {
    get_par("b_mid_Intercept") +
      get_par("b_mid_temp_c") * temp_c_value
  } else if (stage_name == "old_embryos") {
    get_par("b_mid_Intercept") +
      get_par("b_mid_life_stageold_embryos") +
      (get_par("b_mid_temp_c") + get_par("b_mid_temp_c:life_stageold_embryos")) * temp_c_value
  } else if (stage_name == "young_embryos") {
    get_par("b_mid_Intercept") +
      get_par("b_mid_life_stageyoung_embryos") +
      (get_par("b_mid_temp_c") + get_par("b_mid_temp_c:life_stageyoung_embryos")) * temp_c_value
  }
}

get_stage_4pl_draws <- function(stage_name, temp_c_value) {
  data.frame(
    .draw = dr$.draw,
    life_stage = stage_name,
    temp_c = temp_c_value,
    low = if (stage_name == "larvae") {
      plogis(get_par("b_lraw_life_stagelarvae")) * 0.49
    } else if (stage_name == "old_embryos") {
      plogis(get_par("b_lraw_life_stageold_embryos")) * 0.49
    } else {
      plogis(get_par("b_lraw_life_stageyoung_embryos")) * 0.49
    },
    up = if (stage_name == "larvae") {
      0.51 + plogis(get_par("b_uraw_life_stagelarvae")) * 0.49
    } else if (stage_name == "old_embryos") {
      0.51 + plogis(get_par("b_uraw_life_stageold_embryos")) * 0.49
    } else {
      0.51 + plogis(get_par("b_uraw_life_stageyoung_embryos")) * 0.49
    },
    k = if (stage_name == "larvae") {
      exp(get_par("b_logk_life_stagelarvae"))
    } else if (stage_name == "old_embryos") {
      exp(get_par("b_logk_life_stageold_embryos"))
    } else {
      exp(get_par("b_logk_life_stageyoung_embryos"))
    },
    mid = get_mid_draws(stage_name, temp_c_value)
  )
}

stage_levels <- levels(d$life_stage)

# =========================
# 5) TDT curve for each life stage
# =========================
temps_df <- data.frame(
  assay_temp = seq(38, 44, by = 0.01)
) %>%
  mutate(
    temp_c = assay_temp - temp_mean
  )

tdt_draws_avg <- bind_rows(
  lapply(stage_levels, function(stage_name) {
    bind_rows(
      lapply(seq_len(nrow(temps_df)), function(i) {
        out <- get_stage_4pl_draws(
          stage_name = stage_name,
          temp_c_value = temps_df$temp_c[i]
        )
        out$assay_temp <- temps_df$assay_temp[i]
        out
      })
    )
  })
) %>%
  mutate(
    log10_t_h = mid + log((up - 0.5) / (0.5 - low)) / k,
    t50_h = 10^log10_t_h,
    t50_min = 60 * t50_h
  ) %>%
  filter(
    is.finite(t50_min),
    t50_min > 0
  )

pred_df_tdt_avg <- tdt_draws_avg %>%
  group_by(life_stage, assay_temp) %>%
  summarise(
    time_to_50 = median(t50_min, na.rm = TRUE),
    lower = quantile(t50_min, 0.025, na.rm = TRUE),
    upper = quantile(t50_min, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

p_tdt_avg_lin <- ggplot(pred_df_tdt_avg) +
  geom_ribbon(
    aes(x = assay_temp, ymin = lower, ymax = upper, fill = life_stage),
    alpha = 0.20,
    show.legend = FALSE
  ) +
  geom_line(
    aes(x = assay_temp, y = time_to_50, colour = life_stage),
    linewidth = 1.6
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_colour_manual(
    values = life_stage_cols,
    breaks = names(life_stage_labels),
    labels = life_stage_labels
  ) +
  scale_fill_manual(values = life_stage_cols) +
  scale_y_continuous(breaks = ybreaks_left, labels = lab_number) +
  coord_cartesian(ylim = c(0.9, 960)) +
  labs(
    x = "Assay temperature (°C)",
    y = "Time to 50% survival (min)",
    colour = "Life stage"
  ) +
  base_theme +
  theme(
    legend.position = c(0.8, 0.5),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

p_tdt_avg_log <- ggplot(pred_df_tdt_avg) +
  geom_ribbon(
    aes(x = assay_temp, ymin = lower, ymax = upper, fill = life_stage),
    alpha = 0.20,
    show.legend = FALSE
  ) +
  geom_line(
    aes(x = assay_temp, y = time_to_50, colour = life_stage),
    linewidth = 1.6
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_colour_manual(
    values = life_stage_cols,
    breaks = names(life_stage_labels),
    labels = life_stage_labels
  ) +
  scale_fill_manual(values = life_stage_cols) +
  scale_y_log10(breaks = ybreaks_right, labels = lab_number) +
  coord_cartesian(ylim = c(0.9, 960)) +
  labs(
    x = "Assay temperature (°C)",
    y = expression(log[10] * "(Time to 50% survival, min)"),
    colour = "Life stage"
  ) +
  base_theme +
  theme(
    legend.position = "none"
  )

print(p_tdt_avg_lin + p_tdt_avg_log + patchwork::plot_layout(ncol = 2))
