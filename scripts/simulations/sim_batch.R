# Non-independence simulation: scenario 10 of the paper's simulations (Reviewer 2: pseudoreplication and
# nested random effects). Notes: notes/case_studies/2026-11-28-temperature-lack-of-fit.qmd.
# Each day shifts CTmax for every cup it holds by a random amount, drawn anew for each dataset: SD
# 0.53 C of CTmax, which is 0.08 log10 min on the time axis (scen10_A_batch16 doubles it). Schedules as
# in fig_batch_sequence.R. A runs one temperature per day, so a day is that temperature's run. B runs
# every temperature every day (replicate cup r on day r). C runs each temperature on two days (one
# temperature a day, its 5 cups per exposure split 3 + 2). All three test the same 150 cups
# (5 temperatures x 6 durations x 5 cups). Figures: fig_batch_design.R and fig_batch_sequence.R.
# The asymptotes sit at the bounds (0.999 and 0.001, as in scenario 2), so every method is unbiased
# without day effects and the cells isolate what ignoring days does.
#
# Per dataset, three fits:
#   joint      direct-form joint 4PL, CTmax and z as parameters (shapes vary with temperature,
#              as in the paper's simulations)
#   joint_day  the same with a random intercept for day on CTmax, which also estimates the SD of the
#              day shifts. In A a day holds one temperature, so this is the A-design term (one shift
#              per temperature); in B it is the day term (one shift per day, shared by all temperatures);
#              in C it can compare the two days at each temperature
#   two_stage  the classical two-stage estimate with t-intervals
# z and CTmax_1hr are read at 50% survival and scored against the straight TDT line. All fits use the
# package's default priors. The day SD has an exponential(2) prior (mean 0.5 C), close to the true
# 0.53 C, so its coverage alone says little. Each row therefore also stores the SD of the shifts actually
# drawn in that dataset (sd_realised, C of CTmax), to compare with the estimate.
suppressPackageStartupMessages({library(bayesTLS); library(dplyr)})
source(here::here("scripts", "simulations", "sim_functions.R"))

OUT     <- here::here("output", "sim_batch")
N_SIMS  <- 1000         # datasets per scenario, as in the main simulation (2 for a quick pilot)
WORKERS <- 4            # parallel R sessions, 3 chains each (1 = run here, in order)
SEED    <- 30000000L    # dataset seed = SEED + 1000 * scenario number + sim_id; this range is disjoint
                        # from the main simulation's seeds (20,261,514 to 20,290,513)
bounds  <- list(u_0 = 0.999, ell_0 = 0.001)                        # asymptotes at the bounds
truths  <- list(scen10_A_none    = do.call(sim_truth, c(bounds, schedule = "A")),
                scen10_A_batch   = do.call(sim_truth, c(bounds, schedule = "A", batch_sd = 0.08)),
                scen10_B_none    = do.call(sim_truth, c(bounds, schedule = "B")),
                scen10_B_batch   = do.call(sim_truth, c(bounds, schedule = "B", batch_sd = 0.08)),
                scen10_C_none    = do.call(sim_truth, c(bounds, schedule = "C")),
                scen10_C_batch   = do.call(sim_truth, c(bounds, schedule = "C", batch_sd = 0.08)),
                scen10_A_batch16 = do.call(sim_truth, c(bounds, schedule = "A", batch_sd = 0.16)))

# Workers source a frozen copy of the harness, so editing sim_functions.R during the run cannot change
# later datasets. Its md5 is stored in every row.
dir.create(OUT, showWarnings = FALSE)
frozen  <- file.path(OUT, "sim_functions_frozen.R")
invisible(file.copy(here::here("scripts", "simulations", "sim_functions.R"), frozen, overwrite = TRUE))
HARNESS <- unname(tools::md5sum(frozen))

# Pre-flight on seed 1: the days each schedule makes, and different data wherever the truth differs
for (nm in names(truths)) {
  d <- sim_dataset(n_reps = 5, seed = 1L, truth = truths[[nm]])
  cat(sprintf("%-17s z %.2f  CTmax_1hr %.2f | %2d days, up to %d temperatures a day | y hash %s\n", nm,
              truths[[nm]]$z_true, truths[[nm]]$CTmax_1hr_true, length(unique(d$day)),
              max(tapply(d$T, d$day, function(v) length(unique(v)))), substr(rlang::hash(d$y), 1, 10)))
}

# One simulated dataset, start to finish. Resume-safe: a finished dataset (sim_*.rds) is skipped; a
# failed one is saved as fail_*.rds with its error and tried again on the next run.
one_sim <- function(sim_id, scenario) {
  f <- file.path(OUT, "raw", scenario, sprintf("sim_%05d.rds", sim_id))
  if (file.exists(f)) return(invisible(NULL))
  truth   <- truths[[scenario]]
  seed    <- SEED + 1000L * match(scenario, names(truths)) + sim_id
  d       <- sim_dataset(n_reps = 5, seed = seed, truth = truth)
  sd_real <- sd(d$shift[!duplicated(d$day)]) * truth$z_true        # day shifts drawn, in C of CTmax
  row <- tryCatch({
    std  <- standardize_data(d, temp = "T", duration = "t", n_total = "n", n_surv = "y",
                             duration_unit = "minutes")
    fits <- list(
      joint     = fit_4pl(std, ctmax = ~ 1, z = ~ 1, low = ~ temp_c, up = ~ temp_c, k = ~ temp_c,
                          t_ref = 60, seed = seed, chains = 3, iter = 3000, warmup = 1500, cores = 3,
                          refresh = 0, silent = 2, control = list(adapt_delta = 0.99, max_treedepth = 16)),
      joint_day = fit_4pl(std, ctmax = ~ 1 + (1 | day), z = ~ 1, low = ~ temp_c, up = ~ temp_c, k = ~ temp_c,
                          t_ref = 60, seed = seed, chains = 3, iter = 3000, warmup = 1500, cores = 3,
                          refresh = 0, silent = 2, control = list(adapt_delta = 0.99, max_treedepth = 16)))
    est <- bind_rows(lapply(fits, function(w) get_tls_est(tls(w, target_surv = "absolute", t_ref = 60))),
                     .id = "method") |>
      filter(quantity %in% c("z", "CTmax")) |>
      select(method, quantity, estimate = median, lower, upper)
    sd_day <- posterior::summarise_draws(
      posterior::subset_draws(posterior::as_draws_df(get_brmsfit(fits$joint_day)), variable = "^sd_day__", regex = TRUE),
      estimate = median, ~ quantile(.x, c(0.025, 0.975)))
    ts  <- fit_two_stage(d, stage1 = "binomial")
    est <- bind_rows(est,
                     tibble(method = "joint_day", quantity = "sd_day", estimate = sd_day$estimate,
                            lower = sd_day$`2.5%`, upper = sd_day$`97.5%`),
                     tibble(method = "two_stage", quantity = c("z", "CTmax"),
                            estimate = c(ts$z$point, ts$CTmax_1hr$point),
                            lower = c(ts$z$lower_t, ts$CTmax_1hr$lower_t),
                            upper = c(ts$z$upper_t, ts$CTmax_1hr$upper_t)))
    diag <- bind_rows(lapply(fits, diagnose_tdt_fit), .id = "method") |>
      select(method, divergences, rhat_max, all_pass)
    est |>
      left_join(diag, by = "method") |>
      mutate(sim_id = sim_id, scenario = scenario, sd_realised = sd_real, harness = HARNESS,
             truth = case_when(quantity == "z" ~ truth$z_true, quantity == "CTmax" ~ truth$CTmax_1hr_true,
                               quantity == "sd_day" ~ (truth$batch_sd %||% 0) * truth$z_true),
             bias = estimate - truth,
             covered = if_else(quantity == "sd_day" & truth == 0, NA, lower <= truth & truth <= upper),
             width = upper - lower, success = is.finite(estimate))
  }, error = function(e) e)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  if (inherits(row, "error")) {
    save_retry(list(scenario = scenario, sim_id = sim_id, error = conditionMessage(row)),
               file.path(dirname(f), sprintf("fail_%05d.rds", sim_id)))
  } else {
    save_retry(list(row = row), f)
  }
  invisible(NULL)
}

# All datasets still to run, interleaved across scenarios so every worker covers all of them
todo <- bind_rows(lapply(names(truths), function(s)
  tibble(scenario = s, sim_id = setdiff(seq_len(N_SIMS), done_ids(file.path(OUT, "raw", s)))))) |>
  arrange(sim_id)
cat(nrow(todo), "datasets to run\n")
t0 <- Sys.time()
run_parallel(seq_len(nrow(todo)), function(k) one_sim(todo$sim_id[k], todo$scenario[k]), workers = WORKERS,
             exports = c("one_sim", "truths", "todo", "OUT", "SEED", "HARNESS"), src = frozen)
cat(sprintf("finished in %.1f h\n", as.numeric(difftime(Sys.time(), t0, units = "hours"))))

# Summary: bias (mean and median), coverage of the 95% intervals and their width, with Monte Carlo SEs.
# summ_conv repeats it on the datasets where both joint fits had no divergences and R-hat < 1.01, the
# same datasets for all three methods.
fails <- list.files(file.path(OUT, "raw"), pattern = "^fail_", recursive = TRUE)
cat("failed datasets:", length(fails), "\n")
per_sim <- bind_rows(lapply(names(truths), function(s) collect_raw(file.path(OUT, "raw", s))$per_sim))
summ <- summarise_mcse(per_sim) |>
  left_join(per_sim |> group_by(scenario, method, quantity) |>
              summarise(median_bias = median(bias, na.rm = TRUE),
                        share_divergent = mean(divergences > 0),
                        cor_realised = if (first(quantity) == "sd_day" && isTRUE(sd(sd_realised) > 0))
                          cor(estimate, sd_realised) else NA_real_,
                        .groups = "drop"),
            by = c("scenario", "method", "quantity"))
converged <- per_sim |> filter(method != "two_stage") |> group_by(scenario, sim_id) |>
  summarise(ok = all(divergences == 0 & rhat_max < 1.01), .groups = "drop") |> filter(ok)
summ_conv <- summarise_mcse(semi_join(per_sim, converged, by = c("scenario", "sim_id")))
print(as.data.frame(summ))
saveRDS(list(per_sim = per_sim, summary = summ, summary_converged = summ_conv, failed = fails),
        file.path(OUT, "results.rds"))
