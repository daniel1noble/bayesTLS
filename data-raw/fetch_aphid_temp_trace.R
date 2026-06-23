# Re-source the field temperature series for the cereal-aphid heat-injury
# projection (Li et al. 2023). Their Dryad deposit ships only surv.txt; the
# hourly temperature CSVs they obtained from the VisualCrossing API (Wuhan,
# Xinxiang, Beijing; summer 2016) are NOT deposited. We re-source the same
# cities and period from the Open-Meteo historical archive (ERA5 reanalysis;
# free, no API key; https://open-meteo.com/en/docs/historical-weather-api).
#
# ERA5 is a different reanalysis of the same weather, so this trace is comparable
# to -- but not a bit-reproduction of -- Li et al.'s VisualCrossing input. The HI
# projection in the case study uses OUR joint-fit TDT parameters anyway, so it is
# a demonstration of the bayesTLS workflow comparable to their Fig. 2, not an
# exact reproduction. The query is a plain deterministic URL, so anyone can
# re-run this and get the same trace.
#
# Writes inst/extdata/data_temp_trace_aphid_summer2016.csv.
# Run from the package root:  Rscript data-raw/fetch_aphid_temp_trace.R

suppressPackageStartupMessages({
  library(jsonlite); library(dplyr); library(readr); library(tibble)
})

# Wuhan coords from the Li et al. README; Xinxiang/Beijing are city centres.
cities <- tibble::tribble(
  ~city,      ~lat,     ~lon,
  "Wuhan",    30.78,    114.21,
  "Xinxiang", 35.3027,  113.9268,
  "Beijing",  39.9042,  116.4074
)
START <- "2016-05-01"
END   <- "2016-08-31"

fetch_city <- function(city, lat, lon) {
  url <- sprintf(paste0(
    "https://archive-api.open-meteo.com/v1/archive",
    "?latitude=%s&longitude=%s&start_date=%s&end_date=%s",
    "&hourly=temperature_2m&timezone=Asia%%2FShanghai"),
    lat, lon, START, END)
  j <- jsonlite::fromJSON(url)
  tibble::tibble(
    city     = city,
    datetime = j$hourly$time,
    temp_c   = as.numeric(j$hourly$temperature_2m)
  ) |>
    dplyr::mutate(time_h = seq_len(dplyr::n()) - 1L)   # hours from series start
}

trace <- dplyr::bind_rows(
  Map(fetch_city, cities$city, cities$lat, cities$lon)
)
stopifnot(!anyNA(trace$temp_c), nrow(trace) > 0)

readr::write_csv(trace, "inst/extdata/data_temp_trace_aphid_summer2016.csv")
message(sprintf("wrote %d rows (%d cities) to inst/extdata/data_temp_trace_aphid_summer2016.csv",
                nrow(trace), dplyr::n_distinct(trace$city)))
print(trace |> dplyr::group_by(city) |>
        dplyr::summarise(n = dplyr::n(),
                         min_c = round(min(temp_c), 1),
                         max_c = round(max(temp_c), 1),
                         hrs_over_34 = sum(temp_c > 34), .groups = "drop"))
