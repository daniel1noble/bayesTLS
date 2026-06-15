## Build the Oersted et al. (2024) / NicheMapR microclimate trace used for the
## Drosophila suzukii heat-injury panel of manuscript Figure 5.
##
## Source workflow:
##   Oersted et al. (2024), "Thermal limits of survival and reproduction depend
##   on stress duration: a case study of Drosophila suzukii", Zenodo record
##   10821572, file microclimate_injury_accumulation.R.
##
## That script generates a full-year 2018 microclimate at:
##   longitude -1.583131, latitude 48.163910 (Rennes, France)
## using NicheMapR::micro_ncep(..., runshade = 0, minshade = 0, Usrhyt = 0.2),
## and then uses TALOC as the shaded local microclimate temperature.
##
## Reproduce from the repository root:
##   Rscript data-raw/build_orsted_microclimate_trace.R
##
## Output:
##   inst/extdata/orsted_2024/orsted2024_nichemapr_rennes_2018_hourly.csv.gz

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(readr)
  library(NicheMapR)
})

OUT_DIR <- here::here("inst", "extdata", "orsted_2024")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

out <- file.path(OUT_DIR, "orsted2024_nichemapr_rennes_2018_hourly.csv.gz")

loc <- c(-1.583131, 48.163910)

micro <- NicheMapR::micro_ncep(
  loc = loc,
  dstart = "01/01/2018",
  dfinish = "31/12/2018",
  runshade = 0,
  minshade = 0,
  Usrhyt = 0.2
)

metout <- as.data.frame(micro$metout)

trace <- metout |>
  dplyr::mutate(
    hour_index = dplyr::row_number() - 1L,
    datetime_utc = as.POSIXct("2018-01-01 00:00:00", tz = "UTC") +
      hour_index * 3600,
    longitude = loc[1],
    latitude = loc[2]
  ) |>
  dplyr::transmute(
    datetime_utc,
    doy = .data$DOY,
    time_min = .data$TIME,
    hour_index,
    hour_of_day = (.data$TIME / 60) %% 24,
    air_temp_c = .data$TAREF,
    micro_temp_c = .data$TALOC,
    longitude,
    latitude,
    source = "Oersted et al. 2024 Zenodo 10821572; NicheMapR micro_ncep TALOC"
  )

readr::write_csv(trace, out)

message(sprintf(
  "wrote %s (%d rows; micro_temp_c range %.1f to %.1f C)",
  out, nrow(trace), min(trace$micro_temp_c, na.rm = TRUE),
  max(trace$micro_temp_c, na.rm = TRUE)
))
