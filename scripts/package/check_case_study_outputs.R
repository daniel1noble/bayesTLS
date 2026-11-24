#!/usr/bin/env Rscript
# Do the package's derived quantities change between two versions on the SAME
# cached fits? For every fit ms/supplement.qmd caches (the case-study models and
# the tutorial/simulation models) this computes a battery of outputs under a
# BASELINE checkout (a git ref; default 51023b9 = bayesTLS 1.1.0 as it stood
# before the alternative-models work, the version the cached fits were made
# with) and under the WORKING TREE, each in its own R process, and compares
# them for exact equality (identical(), no tolerance).
#
#   Rscript scripts/package/check_case_study_outputs.R             # vs 51023b9
#   Rscript scripts/package/check_case_study_outputs.R <git-ref>
#   Rscript scripts/package/check_case_study_outputs.R --compare-only [<git-ref>]
# Results are kept in output/package_checks/ (ignored by git) so a comparison
# can be re-read without recomputing (--compare-only).
#
# Outputs per fit, each captured with its warnings (or its error message):
#   tls():      z, CTmax and (count data) T_crit, relative and absolute (50%
#               survival) thresholds -- summaries and every posterior draw
#   tdt_parameter_table():  low, up, k and the TDT line per group
#   extract_4pl_pars():     every draw of low, up, k, mid at each temperature
#   derive_z(), derive_temperature_for_duration(), derive_tdt_curve():
#               absolute 50% threshold
#   predict_heat_injury():  relative and absolute dose, constant and varying
#               shape (100 draws, seed 1)
# Each output is computed on the centred scale (temp_mean = 0; the strictest
# bit-level test, since adding ~30 absorbs ulp-level CTmax/T_crit differences)
# and again with temp_mean = 30.5 (keys "<fit>@Tbar30.5"), which exposes
# regressions in how temp_mean is applied (the sign of temp - temp_mean, a
# dropped Tbar +). Neither is the paper's own mean temperature: the raw
# temperatures are not stored in the fits. The baseline is exported with
# `git archive` (read-only).
args <- commandArgs(trailingOnly = TRUE)
root <- here::here()

compute_mode <- length(args) >= 1 && identical(args[1], "--compute")
if (compute_mode) {
  pkg_dir <- args[2]; out_file <- args[3]
  suppressPackageStartupMessages({
    library(brms)
    devtools::load_all(pkg_dir, helpers = FALSE, quiet = TRUE)
  })
  source(file.path(root, "scripts", "package", "supplement_fits.R"))
  run <- function(expr) {                # value + warnings, or the error message
    set.seed(1)                           # same RNG state for every output
    w <- character()
    v <- withCallingHandlers(
      tryCatch(expr, error = function(err)
        structure(conditionMessage(err), class = "run_error")),
      warning = function(x) {
        w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
      })
    list(value = v, warnings = w)
  }
  res <- list()
  for (tbar in c(0, 30.5)) for (e in supplement_cached_fits(root, tbar = tbar)) {
    key <- if (tbar == 0) e$name else paste0(e$name, "@Tbar", tbar)
    if (!identical(e$status, "ok")) { res[[key]] <- list(status = e$status); next }
    wf     <- supplement_workflow(e)
    by     <- if (length(wf$meta$group_vars)) wf$meta$group_vars else NULL
    lethal <- identical(attr(e$data, "tdt_meta")$response_type, "count")
    t_ref  <- wf$meta$t_ref %||% 60
    tc     <- sort(unique(e$data$temp))          # on this pass's scale
    ramp   <- data.frame(time = seq(0, 6, by = 0.1),
                         temp = seq(stats::quantile(e$data$temp, 0.1),
                                    stats::quantile(e$data$temp, 0.9),
                                    length.out = 61))
    keep <- function(x) x[intersect(names(x), c("summary", "draws", "local_summary"))]
    hi <- function(target, shape) {
      run(predict_heat_injury(ramp, wf, target_surv = target, ndraws = 100,
                              seed = 1, shape = shape)$summary)
    }
    res[[key]] <- list(
      status     = "ok",
      tls_rel    = run(keep(tls(wf, by = by, params = c("z", "ctmax", if (lethal) "tcrit"),
                                lethal = lethal, t_ref = t_ref, seed = 1))),
      tls_abs    = run(keep(tls(wf, by = by, params = c("z", "ctmax"),
                                target_surv = 0.5, t_ref = t_ref, seed = 1))),
      par_table  = run(tdt_parameter_table(wf, by = by)),
      pars_draws = run(extract_4pl_pars(wf, by = by)),
      z_abs      = run(keep(derive_z(wf, target_surv = 0.5, by = by))),
      temp_abs   = run(keep(derive_temperature_for_duration(
                     wf, exposure_duration = 60, by = by, target_surv = 0.5,
                     temp_grid = seq(min(tc) - 2, max(tc) + 2, length.out = 400)))),
      curve_abs  = run(keep(derive_tdt_curve(
                     wf, temp_grid = unname(stats::quantile(tc, c(0.1, 0.5, 0.9))),
                     target_surv = 0.5, by = by))),
      hi_rel_const  = hi("relative", "constant"),
      hi_abs_const  = hi(0.5, "constant"),
      hi_rel_vary   = hi("relative", "varying"),
      hi_abs_vary   = hi(0.5, "varying"))
    cat("computed", key, "\n")
  }
  res$.version <- as.character(utils::packageVersion("bayesTLS"))
  saveRDS(res, out_file)
  quit(save = "no")
}

## ---- driver: export the baseline, compute both, compare --------------------
compare_only <- length(args) >= 1 && identical(args[1], "--compare-only")
if (compare_only) args <- args[-1]
ref     <- if (length(args) >= 1) args[1] else "51023b9"
out_dir <- file.path(root, "output", "package_checks")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_base <- file.path(out_dir, sprintf("case_study_outputs_baseline_%s.rds", ref))
out_new  <- file.path(out_dir, "case_study_outputs_working_tree.rds")
if (!compare_only) {
  tmp <- tempfile("bayesTLS-baseline-")
  dir.create(tmp)
  status <- system(sprintf("git -C %s archive %s | tar -x -C %s",
                           shQuote(root), shQuote(ref), shQuote(tmp)))
  if (status != 0) stop("could not export git ref ", ref)
  self <- file.path(root, "scripts", "package", "check_case_study_outputs.R")
  for (job in list(c(tmp, out_base, paste("baseline", ref)),
                   c(root, out_new, "working tree"))) {
    cat("== computing outputs with", job[3], "\n")
    st <- system2("Rscript", c(shQuote(self), "--compute", shQuote(job[1]),
                               shQuote(job[2])))
    if (st != 0) stop("computation failed for ", job[3])
  }
}
base <- readRDS(out_base); new <- readRDS(out_new)
cat(sprintf("\nbaseline %s (bayesTLS %s) vs working tree (bayesTLS %s)\n",
            ref, base$.version, new$.version))
fits <- setdiff(names(base), ".version")
outputs <- c("tls_rel", "tls_abs", "par_table", "pars_draws", "z_abs", "temp_abs",
             "curve_abs", "hi_rel_const", "hi_abs_const", "hi_rel_vary",
             "hi_abs_vary")
# Attributes other than the structural ones (names, class, row.names, dim),
# at the top level and one list level down ("summary$model"), and the values
# with them dropped, to tell "the numbers changed" apart from "an attribute
# was added, removed or changed".
extra_attrs <- function(x) {
  setdiff(names(attributes(x)), c("names", "class", "row.names", "dim", "dimnames"))
}
attr_names <- function(x) {
  a <- extra_attrs(x)
  if (is.list(x) && !is.data.frame(x) && !is.null(names(x)))
    a <- c(a, unlist(lapply(names(x), function(n) paste0(n, "$", extra_attrs(x[[n]])))))
  a
}
drop_extra <- function(x) {
  for (n in extra_attrs(x)) attr(x, n) <- NULL
  if (is.list(x) && !is.data.frame(x)) x[] <- lapply(x, drop_extra)
  x
}
verdict <- function(a, b) {
  if (is.null(a) || is.null(b)) return("missing")
  ea <- inherits(a$value, "run_error"); eb <- inherits(b$value, "run_error")
  if (ea || eb) return(if (ea && eb && identical(unclass(a$value), unclass(b$value)))
                         "same error (nothing compared)" else "ERROR DIFFERS")
  w <- if (identical(a$warnings, b$warnings)) "" else "; warnings differ"
  if (identical(a$value, b$value))
    return(if (nzchar(w)) "same values, warnings differ" else "identical")
  if (identical(drop_extra(a$value), drop_extra(b$value))) {
    added   <- setdiff(attr_names(b$value), attr_names(a$value))
    removed <- setdiff(attr_names(a$value), attr_names(b$value))
    return(sprintf("same values, attributes differ (added: %s; removed: %s%s)%s",
                   paste(added, collapse = ","), paste(removed, collapse = ","),
                   if (!length(added) && !length(removed)) "; values of attributes changed" else "",
                   w))
  }
  paste("VALUES DIFFER:", all.equal(a$value, b$value)[1])
}
tab <- do.call(rbind, lapply(fits, function(f) {
  if (!identical(base[[f]]$status, "ok"))
    return(data.frame(fit = f, output = "-", verdict = base[[f]]$status))
  data.frame(fit = f, output = outputs,
             verdict = vapply(outputs, function(o) verdict(base[[f]][[o]], new[[f]][[o]]), ""))
}))
wide <- stats::reshape(tab, idvar = "fit", timevar = "output", direction = "wide")
names(wide) <- sub("^verdict\\.", "", names(wide))
print(wide, row.names = FALSE, right = FALSE)
attr_only <- startsWith(tab$verdict, "same values, attributes differ (added: model; removed: )")
cat(sprintf(paste0("\n%d comparisons: %d identical, %d identical values whose only ",
                   "change is the added \"model\" attribute, %d other\n"), nrow(tab),
            sum(tab$verdict == "identical"), sum(attr_only),
            sum(tab$verdict != "identical" & !attr_only)))
other <- tab[tab$verdict != "identical", ]
if (nrow(other)) {
  cat("\nNot identical:\n")
  print(other, row.names = FALSE, right = FALSE)
}

## z and CTmax side by side (relative threshold), so the inferences can be read
cat("\nz and CTmax (relative threshold; CTmax on the centred scale), median [95% CrI]\n")
fmt <- function(s, q) {
  r <- s[s$quantity == q, , drop = FALSE]
  if (!nrow(r)) return(character())
  lab <- if (ncol(r) > 0 && any(!names(r) %in% c("quantity", "median", "lower", "upper",
                                                   "t_ref", "target_surv", "mode")))
    apply(r[, setdiff(names(r), c("quantity", "median", "lower", "upper", "t_ref",
                                  "target_surv", "mode")), drop = FALSE], 1,
          paste, collapse = "/") else "all"
  sprintf("%s %s: %.4f [%.4f, %.4f]", q, lab, r$median, r$lower, r$upper)
}
for (f in fits) {
  if (grepl("@Tbar", f, fixed = TRUE)) next     # centred pass shown
  b <- base[[f]]$tls_rel$value; n <- new[[f]]$tls_rel$value
  if (is.null(b) || inherits(b, "run_error") || is.null(n) ||
      inherits(n, "run_error")) next
  rows_b <- c(fmt(b$summary, "z"), fmt(b$summary, "CTmax"))
  rows_n <- c(fmt(n$summary, "z"), fmt(n$summary, "CTmax"))
  cat(sprintf("  %-38s %s\n", f, ""))
  for (i in seq_along(rows_b))
    cat(sprintf("      baseline %-52s | new %s\n", rows_b[i], rows_n[i]))
}
