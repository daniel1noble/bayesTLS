# Build the merged natural-temperature series used in the heat-injury
# demonstration, from the Dryad deposit of Rezende et al. (2020), Science.
#
#   Rezende, Bozinovic, Szilagyi & Santos (2020) Science 369:1242-1245.
#   Data (CC0): https://doi.org/10.5061/dryad.stqjq2c1r  (Dryad version 77627)
#
# The deposit stores air temperature in 13 per-year files in two different
# formats: 1984-1991 are single-station semicolon exports from the Chilean DMC
# (station 330021, Santiago); 2014-2018 are five-site comma exports from the
# INIA agro-meteorological network. This script downloads the public version
# zip, parses both formats, and writes one tidy long file:
#
#   inst/extdata/rezende_2020/rezende2020_chile_hourly_temps.csv.gz
#   columns: datetime, site, temp_c, source
#
# Run from the repository root:  Rscript data-raw/build_rezende_temps.R

suppressPackageStartupMessages({library(dplyr); library(tibble); library(here)})

dest_dir <- here::here("inst", "extdata", "rezende_2020")
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Download + extract the Dryad version zip to a temp dir ----
zip_url <- "https://datadryad.org/api/v2/versions/77627/download"
tmp <- tempfile(fileext = ".zip"); ex <- tempfile()
utils::download.file(zip_url, tmp, mode = "wb", quiet = TRUE)
utils::unzip(tmp, exdir = ex)

site_map <- c(La_Platina_Stgo = "Santiago", Santa_Rosa_Ch = "Chillan",
              Los_Canelos_PM = "Puerto_Montt", Las_Lomas_Val = "Valdivia",
              Curico_Los_Niches = "Curico")

# 1984-1991: CodigoNacional;momento;Ts_Valor (single Santiago DMC station)
parse_early <- function(path) {
  ln <- gsub("\r", "", readLines(path, warn = FALSE)); ln <- ln[nzchar(ln)][-1]
  p  <- do.call(rbind, strsplit(ln, ";", fixed = TRUE))
  tibble(datetime = as.POSIXct(p[, 2], format = "%d-%m-%Y %H:%M:%S", tz = "America/Santiago"),
         site = "Santiago", temp_c = suppressWarnings(as.numeric(p[, 3])),
         source = paste0("DMC-", trimws(p[, 1])))
}
# 2014-2018: row,date,<5 INIA site columns>
parse_recent <- function(path) {
  ln <- gsub("\r", "", readLines(path, warn = FALSE))
  cols <- strsplit(ln[grep("La_Platina", ln)[1]], ",", fixed = TRUE)[[1]]
  site_cols <- cols[nzchar(cols)][-(1:2)]
  dat <- ln[grepl("^[0-9]+,[0-9]{1,2}-[0-9]{1,2}-[0-9]{4} ", ln)]
  p   <- strsplit(dat, ",", fixed = TRUE)
  dt  <- as.POSIXct(vapply(p, `[`, "", 2), format = "%d-%m-%Y %H:%M", tz = "America/Santiago")
  bind_rows(lapply(seq_along(site_cols), function(j)
    tibble(datetime = dt, site = unname(site_map[site_cols[j]]),
           temp_c = suppressWarnings(as.numeric(vapply(p, `[`, "", 2 + j))),
           source = paste0("INIA-", site_cols[j]))))
}

f <- function(y) file.path(ex, sprintf("Temp_%d.csv", y))
all <- bind_rows(
  bind_rows(lapply(1984:1991, function(y) parse_early(f(y)))),
  bind_rows(lapply(2014:2018, function(y) parse_recent(f(y))))
) |>
  filter(!is.na(datetime)) |>
  arrange(site, datetime) |>
  mutate(datetime = format(datetime, "%Y-%m-%d %H:%M:%S"))

# ---- 3. Write one gzip-compressed tidy CSV ----
out <- file.path(dest_dir, "rezende2020_chile_hourly_temps.csv.gz")
con <- gzfile(out, "w"); write.csv(all, con, row.names = FALSE); close(con)
message(sprintf("wrote %s  (%d rows, %.2f MB)", out, nrow(all), file.info(out)$size / 1e6))
