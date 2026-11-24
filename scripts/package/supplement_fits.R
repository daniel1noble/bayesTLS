# Helpers shared by scripts/package/check_cached_fits_stancode.R and
# scripts/package/check_case_study_outputs.R.
#
# supplement_cached_fits() finds every fit ms/supplement.qmd creates with
#   fit_4pl(..., file = file.path(models_dir, "<name>"))
# and returns, per fit, the cached brmsfit, the supplement's own model
# arguments (evaluated) and a data frame rebuilt from the fit's data with the
# metadata standardize_data() would attach (duration unit and response type
# read from the supplement's standardize_data() call for that data object;
# temp_mean = `tbar`, default 0 = the centred scale, because the paper's raw
# temperatures are not stored in the fit). supplement_workflow() turns one entry
# into a bayes_tls workflow with the cached fit attached.
#
# The functions call whichever bayesTLS is loaded (library() or
# devtools::load_all()) at the time they run, so the same code serves a
# baseline checkout and the working tree.

supplement_parse <- function(root) {
  lines  <- readLines(file.path(root, "ms", "supplement.qmd"), warn = FALSE)
  starts <- grep("^```\\{r", lines)
  ends   <- grep("^```\\s*$", lines)
  calls <- list(); assigns <- list()
  for (s in starts) {
    e <- ends[ends > s][1]
    if (is.na(e)) next
    code <- paste(lines[(s + 1):(e - 1)], collapse = "\n")
    ex <- tryCatch(parse(text = code, keep.source = TRUE), error = function(e) NULL)
    if (is.null(ex)) next
    for (x in as.list(ex)) {                    # top-level `sym <- expr`
      if (is.call(x) && identical(x[[1]], as.name("<-")) && is.name(x[[2]]))
        assigns[[as.character(x[[2]])]] <- x[[3]]
    }
    pd <- utils::getParseData(ex)
    for (id in pd$id[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "fit_4pl"]) {
      call_id <- pd$parent[pd$id == pd$parent[pd$id == id]]
      calls[[length(calls) + 1]] <- list(
        line = s + pd$line1[pd$id == id],
        expr = str2lang(utils::getParseText(pd, call_id)))
    }
  }
  list(calls = calls, assigns = assigns)
}

supplement_cached_fits <- function(root,
                                   models_dir = file.path(root, "output", "models"),
                                   tbar = 0) {
  sp <- supplement_parse(root)
  not_model   <- c("chains", "iter", "warmup", "cores", "seed", "control",
                   "file", "file_refit", "refresh", "backend", "init", "silent",
                   "fit")
  to_stancode <- c("sample_prior", "threads", "normalize", "stanvars")
  sd_formals  <- formals(standardize_data)
  out <- list()
  for (cl in sp$calls) {
    a <- as.list(cl$expr)[-1]
    f <- a$file
    if (is.null(f) || !is.call(f) || !identical(f[[2]], as.name("models_dir"))) next
    name  <- f[[3]]
    path  <- file.path(models_dir, paste0(name, ".rds"))
    entry <- list(name = name, line = cl$line, status = "ok")
    if (!file.exists(path)) {
      entry$status <- "no cached file"
      out[[name]] <- entry
      next
    }
    fit <- readRDS(path)
    if (inherits(fit, "bayes_tls")) fit <- fit$fit
    nm  <- names(a)
    if (is.null(nm)) nm <- rep("", length(a))
    data_expr <- if ("data" %in% nm) a[["data"]] else a[[which(nm == "")[1]]]
    data_sym  <- as.character(data_expr)
    # duration unit and response type from the supplement's standardize_data()
    sd_call <- sp$assigns[[data_sym]]
    unit <- sd_formals$duration_unit
    rtype <- "count"
    if (!is.null(sd_call) && is.call(sd_call)) {
      m <- tryCatch(match.call(standardize_data, sd_call), error = function(e) NULL)
      if (!is.null(m$duration_unit)) unit <- eval(m$duration_unit)
      if (!is.null(m$proportion)) rtype <- "proportion"
    }
    d <- fit$data
    d$duration <- 10^d$logd           # model time units, as standardize_data()
    d$temp     <- d$temp_c + tbar     # tbar = 0: the centred scale
    attr(d, "tdt_meta") <- list(temp_mean = tbar, duration_unit = unit,
                                random_effects = NULL, response_type = rtype,
                                response_var = if (rtype == "proportion") "survival")
    env <- new.env(parent = globalenv())
    assign("data", d, envir = env)
    # Evaluate an argument; only when that fails (an object defined elsewhere
    # in the supplement, e.g. `prior = priors_sex`) resolve the missing symbols
    # from the supplement's own assignments, with the call's data object
    # replaced by `data`. Formulas evaluate without lookups.
    resolve <- function(expr) {
      res <- tryCatch(eval(expr, env), error = function(e) e)
      if (!inherits(res, "error")) return(res)
      for (s in all.names(expr, functions = FALSE)) {
        if (!exists(s, envir = env) && !is.null(sp$assigns[[s]])) {
          rhs <- do.call(substitute, list(sp$assigns[[s]],
                                          stats::setNames(list(as.name("data")),
                                                          data_sym)))
          assign(s, eval(rhs, env), envir = env)
        }
      }
      eval(expr, env)
    }
    keep <- setdiff(nm[nzchar(nm)], c("data", not_model))
    args <- tryCatch(stats::setNames(lapply(keep, function(k) resolve(a[[k]])), keep),
                     error = function(e) e)
    if (inherits(args, "error")) {
      entry$status <- paste("UNRESOLVED:", conditionMessage(args))
      out[[name]] <- entry
      next
    }
    entry$fit       <- fit
    entry$data      <- d
    entry$args      <- args[setdiff(names(args), to_stancode)]
    entry$stan_args <- args[intersect(names(args), to_stancode)]
    out[[name]] <- entry
  }
  out
}

# A bayes_tls workflow for one entry: the loaded package's fit_4pl() spec for
# the supplement's call, with the cached brmsfit attached.
supplement_workflow <- function(entry) {
  wf <- do.call(fit_4pl, c(list(entry$data), entry$args, list(fit = FALSE)))
  wf$fit <- entry$fit
  wf
}
