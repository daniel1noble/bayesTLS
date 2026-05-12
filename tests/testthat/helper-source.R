# Source the full R/ function library once at the top of every testthat run.
# Avoids each test file having to repeat the source() boilerplate.

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(posterior)
})

.project_root <- tryCatch(
  here::here(),
  error = function(e) normalizePath(file.path(getwd(), "..", ".."),
                                    winslash = "/", mustWork = FALSE)
)

.lib_files <- c(
  "utils.R", "standardize_data.R", "priors.R", "fit_4pl.R",
  "predict_survival_curves.R", "extract_tdt.R", "tdt_landscape.R",
  "repair.R", "temperature_scenarios.R", "predict_heat_injury.R",
  "plotting.R", "diagnostics.R"
)
for (f in .lib_files) source(file.path(.project_root, "R", f))
