# Test helper: load the bayesTLS package and its hard runtime deps.
#
# `devtools::test()` and `R CMD check` load the package via
# pkgload::load_all() before running any tests, which exposes both exported
# and internal functions to the test namespace. This file handles the
# plain `testthat::test_dir("tests/testthat")` case by falling back to
# pkgload::load_all() when no namespace is loaded yet.

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(posterior)
})

if (!"bayesTLS" %in% loadedNamespaces()) {
  if (requireNamespace("pkgload", quietly = TRUE) &&
      requireNamespace("here", quietly = TRUE)) {
    pkgload::load_all(here::here(), quiet = TRUE)
  } else {
    library(bayesTLS)
  }
}
