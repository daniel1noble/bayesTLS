#!/usr/bin/env Rscript
# Phase-0 gate for the alternative-models work: for every cached bayes_tls
# workflow in output/models/, regenerate the Stan code from the STORED formula,
# priors and data with the CURRENTLY installed brms, and compare it to the Stan
# code the fit was actually compiled from. A mismatch means brms (or the stored
# formula/prior objects) would make `file_refit = "on_change"` refit that model
# on the next render. Read-only; prints one row per cached fit.
#
#   Rscript scripts/package/check_cached_fits_stancode.R
suppressPackageStartupMessages({ library(bayesTLS); library(brms) })
files <- list.files(here::here("output", "models"), pattern = "\\.rds$", full.names = TRUE)
rows <- lapply(files, function(f) {
  x <- tryCatch(readRDS(f), error = function(e) NULL)
  # fit_4pl(file = ) hands `file` to brms::brm(), so the cache holds the
  # brmsfit itself; a bayes_tls workflow wraps one in $fit.
  fit <- if (inherits(x, "brmsfit")) x else if (inherits(x, "bayes_tls")) x$fit else NULL
  if (is.null(fit))
    return(data.frame(file = basename(f), status = "not a brms fit"))
  stored <- tryCatch(brms::stancode(fit), error = function(e) NA_character_)
  now    <- tryCatch(suppressWarnings(
                       brms::make_stancode(fit$formula, data = fit$data,
                                           prior = fit$prior, data2 = fit$data2)),
                     error = function(e) paste("ERROR:", conditionMessage(e)))
  data.frame(file = basename(f),
             status = if (identical(stored, now)) "identical" else
               if (startsWith(now, "ERROR")) now else "DIFFERS")
})
res <- do.call(rbind, rows)
print(res, row.names = FALSE, right = FALSE)
cat(sprintf("\n%d fitted workflows: %d identical, %d differ, %d errors\n",
            sum(res$status != "not a brms fit"),
            sum(res$status == "identical"), sum(res$status == "DIFFERS"),
            sum(startsWith(res$status, "ERROR"))))
