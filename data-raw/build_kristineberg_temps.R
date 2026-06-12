## Build the Kristineberg (Gullmar fjord) hourly sea-temperature dataset used by
## the brown shrimp (Crangon crangon) panel of manuscript Figure 5.
##
## Source: Kristineberg Marine Research Station automatic weather station,
##   https://www.weather.mi.gu.se/kristineberg/en/data.shtml
##   (Sven Lovén Centre / University of Gothenburg, Gullmar fjord, Sweden).
##   The `result.php` form returns 5-minute records as CSV. We pull the full
##   12-variable field set (matching P. Pottier's earlier processing) for every
##   complete year on record, then keep the continuous `sea_temp` (and
##   `air_temp`) sensor and aggregate to hourly means.
##
## DATA-USE CAVEAT (from the station's data page): these measurements lack
## formal quality assurance and may contain errors; the station asks to be
## contacted before any use in a publication. Confirm with the station before
## the figure goes into the submitted manuscript.
##
## Reproduce:  Rscript data-raw/build_kristineberg_temps.R
## Output:     inst/extdata/kristineberg/kristineberg_sea_temp_hourly.csv.gz

suppressPackageStartupMessages({
  library(dplyr); library(lubridate); library(readr); library(here)
})

YEARS    <- 2008:2024
RAW_DIR  <- here("data-raw", "kristineberg")
OUT_DIR  <- here("inst", "extdata", "kristineberg")
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## All years are stored stacked in ONE gzipped 5-minute archive (the per-year
## downloads are merged and removed). Re-download only if that archive is gone.
MERGED <- file.path(RAW_DIR, "kristineberg_field_temp_5min_all.csv.gz")

## ---- 1. Download (one request per year, all 12 fields) then merge ---------
## The form POSTs date0/date1 + a repeated field[] list + format to result.php.
if (!file.exists(MERGED)) {
  url    <- "https://www.weather.mi.gu.se/kristineberg/en/result.php"
  fields <- c("localtime", "SBE19-01-salinity", "SBE19-01-temperature", "AirTemp",
              "SeaLevel", "SeaTemp", "Salinity32m", "Temp32m", "Salinity1m",
              "Salinity5m", "Temp5m", "Temp1m")
  field_q <- paste0("field%5B%5D=", fields, collapse = "&")
  tmp_dir <- file.path(RAW_DIR, "raw"); dir.create(tmp_dir, showWarnings = FALSE)
  for (yr in YEARS) {
    f <- file.path(tmp_dir, sprintf("field_temp_%d.csv", yr))
    body <- sprintf("date0=%d-01-01&time0=00:00&date1=%d-12-31&time1=23:59&%s&format=text-csv",
                    yr, yr, field_q)
    system2("curl", c("-s", shQuote(url), "--data", shQuote(body),
                      "-o", shQuote(f)))
    message(sprintf("downloaded %d (%s)", yr, format(file.info(f)$size)))
  }
  # Merge the yearly files into one archive (single header) and drop the parts.
  yfiles <- file.path(tmp_dir, sprintf("field_temp_%d.csv", YEARS))
  hdr    <- readLines(yfiles[1], n = 1)
  con    <- gzfile(MERGED, "w"); writeLines(hdr, con)
  for (f in yfiles) writeLines(readLines(f)[-1], con)
  close(con); unlink(tmp_dir, recursive = TRUE)
  message("merged ", length(yfiles), " years -> ", basename(MERGED))
}

## ---- 2. Read + clean (P. Pottier's column mapping) ------------------------
field_temp <- as.data.frame(readr::read_csv(MERGED, show_col_types = FALSE,
                                            na = "no data"))
field_temp <- field_temp[, seq_len(13)]                  # drop trailing-comma col
names(field_temp) <- c(
  "timestamp_utc", "local_time_cet", "salinity_1m_sjoboden", "temp_1m_sjoboden",
  "air_temp", "sea_level_mm", "sea_temp", "deep_water_salinity_lab",
  "deep_water_temp_lab", "salinity_1m_ctd", "surface_water_salinity_lab",
  "surface_water_temp", "temp_1m_ctd")
names(field_temp)[is.na(names(field_temp))] <- "tmp"

field_temp$timestamp_utc      <- as.POSIXct(field_temp$timestamp_utc, tz = "UTC")
field_temp$sea_temp           <- as.numeric(field_temp$sea_temp)
field_temp$air_temp           <- as.numeric(field_temp$air_temp)
field_temp$surface_water_temp <- as.numeric(field_temp$surface_water_temp)

## Drop the lab-sensor error spikes (P. Pottier's filter) and obviously bad
## continuous-sensor values (physical range for this fjord surface).
field_temp <- field_temp |>
  filter(is.na(surface_water_temp) | surface_water_temp < 100) |>
  filter(!is.na(timestamp_utc), !is.na(sea_temp),
         sea_temp > -5, sea_temp < 40)

## ---- 3. Aggregate the 5-min record to hourly means ------------------------
hourly <- field_temp |>
  mutate(datetime_utc = floor_date(timestamp_utc, "hour")) |>
  group_by(datetime_utc) |>
  summarise(sea_temp_c = mean(sea_temp, na.rm = TRUE),
            air_temp_c = mean(air_temp, na.rm = TRUE), .groups = "drop") |>
  mutate(year = year(datetime_utc)) |>
  arrange(datetime_utc) |>
  select(datetime_utc, year, sea_temp_c, air_temp_c)

out <- file.path(OUT_DIR, "kristineberg_sea_temp_hourly.csv.gz")
write_csv(hourly, out)
message(sprintf("wrote %s  (%d hourly rows, %d-%d)",
                out, nrow(hourly), min(hourly$year), max(hourly$year)))

## ---- 4. Per-year summer warmth, to choose the heatwave-scenario year ------
summer <- hourly |>
  filter(month(datetime_utc) %in% 6:8) |>
  group_by(year) |>
  summarise(jja_mean = round(mean(sea_temp_c), 2),
            jja_p95  = round(quantile(sea_temp_c, 0.95), 2),
            jja_max  = round(max(sea_temp_c), 2),
            n_hours  = n(), .groups = "drop") |>
  arrange(desc(jja_mean))
print(as.data.frame(summer))
