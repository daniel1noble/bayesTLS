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
supp  <- file.path(root, "ms", "supplement.qmd")
lines <- readLines(supp, warn = FALSE)
starts <- grep("^```\\{r", lines); ends <- grep("^```\\s*$", lines)
calls <- list(); assigns <- list()
for (s in starts) {
  e <- ends[ends > s][1]; if (is.na(e)) next
  code <- paste(lines[(s + 1):(e - 1)], collapse = "\n")
  ex <- tryCatch(parse(text = code, keep.source = TRUE), error = function(e) NULL)
  if (is.null(ex)) next
  for (x in as.list(ex)) {                      # top-level `sym <- expr`
    if (is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]]))
      assigns[[as.character(x[[2]])]] <- x[[3]]
  }
  pd <- utils::getParseData(ex)
  for (id in pd$id[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "fit_4pl"]) {
    call_id <- pd$parent[pd$id == pd$parent[pd$id == id]]
    calls[[length(calls) + 1]] <- list(line = s + pd$line1[pd$id == id],
                                       expr = str2lang(utils::getParseText(pd, call_id)))
  }
}
not_model <- c("chains", "iter", "warmup", "cores", "seed", "control", "file",
               "file_refit", "refresh", "backend", "init", "silent", "fit")
to_stancode <- c("sample_prior", "threads", "normalize", "stanvars")
std_formals <- formals(standardize_data)
pass2 <- list()
for (cl in calls) {
  a <- as.list(cl$expr)[-1]
  f <- a$file
  if (is.null(f) || !is.call(f) || !identical(f[[2]], as.name("models_dir"))) next
  name <- f[[3]]
  path <- file.path(models_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    pass2[[length(pass2) + 1]] <- data.frame(fit = name, line = cl$line, status = "no cached file")
    next
  }
  fit <- readRDS(path); if (inherits(fit, "bayes_tls")) fit <- fit$fit
  nm  <- names(a); if (is.null(nm)) nm <- rep("", length(a))
  data_expr <- if ("data" %in% nm) a[["data"]] else a[[which(nm == "")[1]]]
  data_sym  <- as.character(data_expr)
  # duration unit and response type from the supplement's standardize_data()
  sd_call <- assigns[[data_sym]]
  unit <- std_formals$duration_unit; rtype <- "count"
  if (!is.null(sd_call) && is.call(sd_call)) {
    m <- tryCatch(match.call(standardize_data, sd_call), error = function(e) NULL)
    if (!is.null(m$duration_unit)) unit <- eval(m$duration_unit)
    if (!is.null(m$proportion)) rtype <- "proportion"
  }
  d <- fit$data
  attr(d, "tdt_meta") <- list(temp_mean = 0, duration_unit = unit,
                              random_effects = NULL, response_type = rtype,
                              response_var = if (rtype == "proportion") "survival")
  env <- new.env(parent = globalenv()); assign("data", d, envir = env)
  # Evaluate an argument; only when that fails (an object defined elsewhere in
  # the supplement, e.g. `prior = priors_sex`) resolve the missing symbols from
  # the supplement's own assignments. Formulas evaluate without lookups.
  resolve <- function(expr) {
    out <- tryCatch(eval(expr, env), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    for (s in all.names(expr, functions = FALSE)) {
      if (!exists(s, envir = env) && !is.null(assigns[[s]])) {
        rhs <- do.call(substitute, list(assigns[[s]],
                                        stats::setNames(list(as.name("data")), data_sym)))
        assign(s, eval(rhs, env), envir = env)
      }
    }
    eval(expr, env)
  }
  keep <- setdiff(nm[nzchar(nm)], c("data", not_model))
  args <- tryCatch(stats::setNames(lapply(keep, function(k) resolve(a[[k]])), keep),
                   error = function(e) e)
  if (inherits(args, "error")) {
    pass2[[length(pass2) + 1]] <- data.frame(fit = name, line = cl$line,
      status = paste("UNRESOLVED:", conditionMessage(args)))
    next
  }
  sc_args <- args[intersect(names(args), to_stancode)]
  wf <- do.call(fit_4pl, c(list(d), args[setdiff(names(args), to_stancode)],
                           list(fit = FALSE)))
  now <- tryCatch(suppressWarnings(do.call(brms::make_stancode,
           c(list(wf$formula, data = fit$data, prior = wf$prior, data2 = fit$data2),
             sc_args))), error = function(e) paste("ERROR:", conditionMessage(e)))
  old <- brms::stancode(fit)
  pass2[[length(pass2) + 1]] <- data.frame(fit = name, line = cl$line,
    status = if (identical(old, now)) "identical" else
      if (startsWith(now, "ERROR")) now else paste("DIFFERS,", first_diff(old, now)))
}
pass2 <- do.call(rbind, pass2)
cat("== Pass 2: supplement fit_4pl() calls rebuilt with the working-tree package\n")
print(pass2, row.names = FALSE, right = FALSE)
cat(sprintf("\nPass 2: %d supplement fits: %d identical, %d differ, %d other\n",
            nrow(pass2), sum(pass2$status == "identical"),
            sum(startsWith(pass2$status, "DIFFERS")),
            sum(!pass2$status %in% "identical" & !startsWith(pass2$status, "DIFFERS"))))
