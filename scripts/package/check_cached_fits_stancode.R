#!/usr/bin/env Rscript
# Two read-only gates for the cached brms fits in output/models/.
#
# Pass 1 -- brms stability. Regenerate each cached fit's Stan code from its own
#   STORED formula, priors and data with the installed brms and compare it with
#   the code the fit was compiled from. A difference means a brms update alone
#   would make file_refit = "on_change" refit that model.
#
# Pass 2 -- builder regression. Rebuild every fit the supplement creates with
#   fit_4pl(..., file = file.path(models_dir, "<name>")) from the supplement's
#   OWN call (parsed from ms/supplement.qmd, so this list cannot drift from the
#   document), using the WORKING-TREE package, and compare the Stan code brms
#   would generate with the cached fit's. A difference means the package now
#   builds a different model for that call, and the next supplement render
#   would refit it. Data come from the cached fit (fit$data); the duration unit
#   and response type come from the supplement's standardize_data() call for
#   the data object; symbols such as `prior = priors_sex` are resolved from the
#   supplement's own assignment, with its data argument replaced.
#
#   Rscript scripts/package/check_cached_fits_stancode.R
suppressPackageStartupMessages(library(brms))
root       <- here::here()
models_dir <- file.path(root, "output", "models")
first_diff <- function(a, b) {
  a <- strsplit(a, "\n")[[1]]; b <- strsplit(b, "\n")[[1]]
  i <- which(a[seq_len(min(length(a), length(b)))] != b[seq_len(min(length(a), length(b)))])[1]
  if (is.na(i)) sprintf("lengths %d vs %d lines", length(a), length(b))
  else sprintf("line %d: cached `%s` | now `%s`", i, trimws(a[i]), trimws(b[i]))
}

## ---- Pass 1: brms stability -------------------------------------------------
files <- list.files(models_dir, pattern = "\\.rds$", full.names = TRUE)
pass1 <- do.call(rbind, lapply(files, function(f) {
  x <- tryCatch(readRDS(f), error = function(e) NULL)
  # fit_4pl(file = ) hands `file` to brms::brm(), so the cache holds the
  # brmsfit itself; a bayes_tls workflow wraps one in $fit.
  fit <- if (inherits(x, "brmsfit")) x else if (inherits(x, "bayes_tls")) x$fit else NULL
  if (is.null(fit)) return(data.frame(file = basename(f), status = "not a brms fit"))
  stored <- tryCatch(brms::stancode(fit), error = function(e) NA_character_)
  now    <- tryCatch(suppressWarnings(
                       brms::make_stancode(fit$formula, data = fit$data,
                                           prior = fit$prior, data2 = fit$data2)),
                     error = function(e) paste("ERROR:", conditionMessage(e)))
  data.frame(file = basename(f),
             status = if (identical(stored, now)) "identical" else
               if (startsWith(now, "ERROR")) now else "DIFFERS")
}))
cat("== Pass 1: stored formula/priors/data -> Stan code with brms",
    as.character(utils::packageVersion("brms")), "\n")
print(pass1, row.names = FALSE, right = FALSE)
n1 <- sum(pass1$status != "not a brms fit")
cat(sprintf("\nPass 1: %d fitted workflows: %d identical, %d differ, %d errors\n\n",
            n1, sum(pass1$status == "identical"), sum(pass1$status == "DIFFERS"),
            sum(startsWith(pass1$status, "ERROR"))))

## ---- Pass 2: builder regression (working-tree package) ----------------------
suppressPackageStartupMessages(devtools::load_all(root, helpers = FALSE, quiet = TRUE))
source(file.path(root, "scripts", "package", "supplement_fits.R"))
pass2 <- lapply(supplement_cached_fits(root, models_dir), function(e) {
  if (!identical(e$status, "ok"))
    return(data.frame(fit = e$name, line = e$line, status = e$status))
  wf  <- supplement_workflow(e)
  now <- tryCatch(suppressWarnings(do.call(brms::make_stancode,
           c(list(wf$formula, data = e$fit$data, prior = wf$prior,
                  data2 = e$fit$data2), e$stan_args))),
           error = function(err) paste("ERROR:", conditionMessage(err)))
  old <- brms::stancode(e$fit)
  data.frame(fit = e$name, line = e$line,
             status = if (identical(old, now)) "identical" else
               if (startsWith(now, "ERROR")) now else
                 paste("DIFFERS,", first_diff(old, now)))
})
pass2 <- do.call(rbind, pass2)
cat("== Pass 2: supplement fit_4pl() calls rebuilt with the working-tree package\n")
print(pass2, row.names = FALSE, right = FALSE)
cat(sprintf("\nPass 2: %d supplement fits: %d identical, %d differ, %d other\n",
            nrow(pass2), sum(pass2$status == "identical"),
            sum(startsWith(pass2$status, "DIFFERS")),
            sum(!pass2$status %in% "identical" & !startsWith(pass2$status, "DIFFERS"))))
