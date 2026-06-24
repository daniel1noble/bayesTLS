# Build the analysis-ready case-study datasets shipped with bayesTLS from the
# raw CSVs in inst/extdata/. The CSVs are the canonical raw source; this script
# documents the (light) cleaning that turns each into the model-ready frame, and
# writes data/<name>.rda via usethis::use_data().
#
# Run from the package root:  Rscript data-raw/make_datasets.R
#
# Datasets created: shrimp_lethal, shrimp_sublethal, zebrafish_lethal,
# snowgum_psii. Documented in R/data.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Resolve a raw file: prefer the in-tree inst/extdata (when running pre-install),
# fall back to the installed package location.
ext <- function(f) {
  p <- file.path("inst", "extdata", f)
  if (file.exists(p)) p else system.file("extdata", f, package = "bayesTLS")
}

## 1. Brown shrimp — lethal TDT -----------------------------------------------
# Keep the columns standardize_data() and the case study use. In the raw CSV
# `Mortality_after_trial` is the death PROPORTION (deaths / N_individuals_after_trial,
# in [0, 1]) — it is consumed downstream by standardize_data(mortality = ...),
# which derives n_surv = round((1 - mortality) * n_total). It must therefore stay
# numeric: a previous as.integer() here floored every proportion < 1 to 0, which
# collapsed the shipped death counts to ~all-zero (fixed 2026-06-16).
shrimp_lethal <- read_csv(ext("data_lethal_TDT_brown_shrimp.csv"),
                          show_col_types = FALSE) |>
  dplyr::transmute(
    Date,
    Tank,
    Temperature_assay        = as.numeric(Temperature_assay),
    Duration_exposure_hours  = as.numeric(Duration_exposure_hours),
    N_individuals_after_trial = as.integer(N_individuals_after_trial),
    Mortality_after_trial    = as.numeric(Mortality_after_trial)
  ) |>
  as.data.frame()

## 2. Brown shrimp — sublethal time-to-knockdown -----------------------------
# Parse the clock-time strings into elapsed minutes to knockdown; drop excluded
# rows. One row per cup.
parse_time_min <- function(s) {
  t <- lubridate::parse_date_time(s, orders = c("I:M:S p", "H:M:S"),
                                  quiet = TRUE)
  lubridate::hour(t) * 60 + lubridate::minute(t) + lubridate::second(t) / 60
}
shrimp_sublethal <- read_csv(ext("data_sublethal_TDT_brown_shrimp.csv"),
                             show_col_types = FALSE) |>
  dplyr::filter(Exclude == "no") |>
  dplyr::mutate(
    assay_temp      = as.numeric(Assay_temperature),
    date_experiment = as.character(lubridate::dmy(Date)),
    tank_ID         = as.character(Tank),
    cup_ID          = paste0(Trial_ID, "_", Sample),
    time_to_event   = parse_time_min(Time_to_death) -
                      parse_time_min(Starting_time)   # minutes
  ) |>
  dplyr::filter(is.finite(time_to_event), time_to_event > 0) |>
  dplyr::select(assay_temp, time_to_event, date_experiment, tank_ID, cup_ID) |>
  as.data.frame()

## 3. Zebrafish — lethal TDT across life stages (counts) ---------------------
# Aggregate the per-day morning/afternoon mortality columns into one death count
# per trial, derive survivors, and keep one row per assay trial.
zf_raw <- read_csv(ext("data_lethal_TDT_zebrafish.csv"), show_col_types = FALSE)
mort_cols <- grep("^Mortality_day_\\d+_(morning|afternoon)$",
                  names(zf_raw), value = TRUE)
zebrafish_lethal <- zf_raw |>
  dplyr::filter(Exclude == "no") |>
  dplyr::mutate(
    n_total    = as.integer(N_total),
    n_dead     = as.integer(rowSums(
      dplyr::across(dplyr::all_of(mort_cols), as.numeric), na.rm = TRUE)),
    n_surv     = pmax(n_total - n_dead, 0L),
    assay_temp = as.numeric(Temperature_assay),
    duration_h = as.numeric(Duration_exposure_hours),
    life_stage = factor(Life_stage,
                        levels = c("young_embryos", "old_embryos", "larvae")),
    Date_experiment = as.character(Date_experiment)
  ) |>
  dplyr::filter(is.finite(duration_h), duration_h > 0,
                is.finite(assay_temp), n_total > 0) |>
  dplyr::select(assay_temp, duration_h, n_total, n_surv, n_dead,
                life_stage, Date_experiment) |>
  as.data.frame()

## 4. Snow gum leaf — PSII functional impairment (continuous proportion) -----
# Arnold et al. (2026) preprint open data (CC BY-NC 4.0; bioRxiv
# 10.64898/2026.04.09.717599), Experiment 1 (light vs dark recovery), snow gum
# (Eucalyptus pauciflora) slice. ~1 cm^2 leaf sections were heat-treated in a
# water bath under sub-saturating light (Temp 30-56 C x Time 5-120 min), then
# left for 90 min in moderate light OR in darkness post-heat (`Recovery_cond`),
# before a final Fv/Fm 16-24 h later. Response is retained PSII function, the
# ratio of post- to pre-heat Fv/Fm (a proportion in [0, 1]); recomputed here from
# the raw final/initial Fv/Fm so it is exact rather than the rounded raw column.
# `Plant_rep` indexes the 6 replicate mature trees; `Measurement_Day` the two
# assay days. Light cleaning only: drop the constant `Species` column, recompute
# `fvfm_prop`, type the grouping factors; no rows dropped.
snowgum_psii <- read_csv(ext("data_function_PSII_TDT_snowgum.csv"),
                         show_col_types = FALSE) |>
  dplyr::transmute(
    Temp            = as.numeric(Temp),        # assay temperature (degrees C)
    Time            = as.numeric(Time),        # exposure duration (minutes)
    recovery        = factor(Recovery_cond,
                             levels = c("Dark", "Light")),  # post-heat light
    plant           = factor(Plant_rep,
                             levels = sort(unique(Plant_rep))),  # 6 trees
    meas_day        = factor(Measurement_Day,
                             levels = sort(unique(Measurement_Day))),  # 2 days
    initial_fvfm    = as.numeric(initial_fvfm),
    final_fvfm      = as.numeric(final_fvfm),
    fvfm_prop       = as.numeric(final_fvfm) / as.numeric(initial_fvfm)
  ) |>
  as.data.frame()

## 5. Drosophila suzukii — multi-trait TDT, per individual -------------------
# Ørsted et al. (2024) open data (CC BY 4.0; Zenodo 10.5281/zenodo.10602268),
# long-format "all_data_long_R3.csv": one row per fly, three thermal-tolerance
# endpoints in a single frame — mortality (`dead`, aggregate to counts), heat
# coma (`t_coma`, a time-to-event; NA where no coma was recorded), and
# reproductive productivity (`prod`, offspring/female/day). Light cleaning only:
# fix column types and order; no rows dropped.
dsuzukii <- read_csv(ext("data_multitrait_TDT_drosophila_suzukii.csv"),
                     show_col_types = FALSE) |>
  dplyr::transmute(
    id     = as.character(id),
    temp   = as.numeric(temp),       # assay temperature (degrees C)
    lvl    = as.numeric(lvl),        # exposure as % of estimated median t_coma
    time   = as.numeric(time),       # exposure duration (minutes)
    sex    = factor(sex, levels = c("F", "M")),
    rep    = as.integer(rep),        # replicate vial within a temp x lvl x sex cell
    prod   = as.numeric(prod),       # offspring / female / day
    dead   = as.integer(dead),       # 1 = died, 0 = survived
    t_coma = as.numeric(t_coma)      # time to heat coma (minutes); NA if none recorded
  ) |>
  as.data.frame()

## 6. Zebrafish larvae — lethal TDT across an oxygen gradient (counts) --------
# Saruhashi et al. (2026) open data (CC BY 4.0; Zenodo 10.5281/zenodo.20075355),
# "Upper thermal limit" sheet. Survival of diploid/triploid larvae assayed at
# 26/38/39/40 C for 3.8-240 min under three oxygen treatments. The categorical
# moderator is `Treatment` (hypoxia/normoxia/hyperoxia); `O` is the nominal % air
# saturation it targets (25/100/225) and `oxygen` is the MEASURED saturation, both
# carried through for reference. `T` is the target assay temperature, `temperature`
# the measured one. Light cleaning only: relabel ploidy, type columns, no rows
# dropped. Full design (all ploidy x oxygen x temp incl. the 26 C control).
zebrafish_o2 <- read_csv(ext("data_lethal_TDT_zebrafish_oxygen.csv"),
                         show_col_types = FALSE) |>
  dplyr::rename(o2_measured = oxygen) |>      # raw measured % air saturation
  dplyr::transmute(
    cohort        = as.character(Cohort),
    ploidy        = factor(Ploidy, levels = c(2, 3),
                           labels = c("diploid", "triploid")),
    oxygen        = factor(Treatment,
                           levels = c("hypoxia", "normoxia", "hyperoxia")),
    o2_nominal    = as.integer(O),            # nominal % air saturation (25/100/225)
    o2_measured   = as.numeric(o2_measured),  # measured % air saturation
    temp          = as.numeric(T),            # target assay temperature (C)
    temp_measured = as.numeric(temperature),  # measured assay temperature (C)
    duration_min  = as.numeric(time),         # exposure duration (minutes)
    n_total       = as.integer(total),
    n_surv        = as.integer(survival)
  ) |>
  as.data.frame()

## 7. Cereal aphids — lethal TDT, three species x three ages (counts) ---------
# Li et al. (2023) open data (CC0; Dryad 10.5061/dryad.mcvdnck4j), surv.txt.
# Survival of three aphid species at three ages across a heat branch (34-40 C)
# and a cold branch (-11 to -3 C); `branch` flags which. Species codes are
# relabelled (MD/SA/RP; the README's "SD" is a typo for the data's SA). Light
# cleaning only: relabel, type, derive `branch`; no rows dropped. Full design.
aphid_tdt <- read_csv(ext("data_lethal_TDT_aphid.csv"), show_col_types = FALSE) |>
  dplyr::transmute(
    species      = factor(dplyr::recode(spec, MD = "M_dirhodum",
                                        SA = "S_avenae", RP = "R_padi"),
                          levels = c("M_dirhodum", "S_avenae", "R_padi")),
    age          = factor(age, levels = c(2, 6, 12)),
    branch       = factor(ifelse(temp >= 0, "heat", "cold"),
                          levels = c("heat", "cold")),
    temp         = as.numeric(temp),          # assay temperature (C)
    duration_min = as.numeric(dur),           # exposure duration (minutes)
    n_total      = as.integer(total),
    n_surv       = as.integer(surv)
  ) |>
  as.data.frame()

usethis::use_data(shrimp_lethal, shrimp_sublethal, zebrafish_lethal,
                  snowgum_psii, dsuzukii, zebrafish_o2, aphid_tdt,
                  overwrite = TRUE)
