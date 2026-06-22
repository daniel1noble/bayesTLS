# Posterior survival predictions on a temperature × duration grid.
# Builds the prediction grid that brms expects, draws population-level
# (random-effects-marginalised) survival probabilities, and summarises them
# with median + 95% credible band per (temp, duration).

#' Build a prediction grid for a fitted TDT workflow
#'
#' Constructs a `newdata` tibble covering all combinations of `temps` and
#' `durations`, adds the derived `logd` and `temp_c` columns the model needs,
#' and fills random-effect grouping columns with their first observed level
#' (these are marginalised out at prediction time via `re_formula = NA`).
#'
#' @param workflow  A `bayes_tls` object.
#' @param temps     Numeric vector of assay temperatures (°C).
#' @param durations Numeric vector of exposure durations (same unit as the
#'                  training data, typically hours).
#' @param n_total   Integer trials count for count families. `NULL` (default)
#'                  uses the median of the training data. Ignored for a
#'                  continuous-proportion (Beta) fit, which has no denominator.
#' @param by        Optional moderator column(s) to cross with the temp x
#'                  duration grid (one block of rows per group, tagged in
#'                  `.grp`). `NULL` (default) builds the ungrouped grid.
#' @return A tibble with one row per `(temp, duration)` pair (per group when
#'         `by` is supplied).
#' @keywords internal
new_tdt_grid <- function(workflow,
                         temps,
                         durations,
                         n_total = NULL,
                         by      = NULL) {
  data <- workflow$data
  meta <- workflow$meta

  # n_total is only meaningful for count families (binomial / beta_binomial). A
  # Beta (continuous-proportion) fit has no denominator, and posterior_linpred()
  # does not need one, so the grid omits it. Detect via metadata, falling back
  # to the presence of an n_total column for hand-built workflows.
  response_type <- meta$response_type %||%
    (if ("n_total" %in% names(data)) "count" else "proportion")
  is_count <- !identical(response_type, "proportion")

  td <- expand.grid(temp     = temps,
                    duration = durations,
                    KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)
  # Grouped fits: cross the moderator level combinations with temp x duration so
  # posterior_linpred() has the moderator column it needs (the single-condition
  # grid is unchanged when by = NULL).
  if (is.null(by)) {
    nd <- td
  } else {
    lev <- unique(data[, by, drop = FALSE])
    nd  <- do.call(rbind, lapply(seq_len(nrow(lev)), function(i) {
      g <- td
      for (v in by) g[[v]] <- rep(lev[[v]][i], nrow(td))
      g
    }))
    rownames(nd) <- NULL
  }

  if (is_count) {
    if (is.null(n_total) && "n_total" %in% names(data)) {
      n_total <- as.integer(round(stats::median(data$n_total, na.rm = TRUE)))
    }
    if (!is.null(n_total)) nd$n_total <- n_total
  }

  nd$logd    <- log10(nd$duration)
  nd$temp_c  <- nd$temp - meta$temp_mean

  for (re_var in tdt_random_effect_variables(meta$random_effects)) {
    if (re_var %in% names(data) && !(re_var %in% names(nd))) {
      if (is.factor(data[[re_var]])) {
        nd[[re_var]] <- factor(levels(data[[re_var]])[1],
                               levels = levels(data[[re_var]]))
      } else {
        nd[[re_var]] <- data[[re_var]][1]
      }
    }
  }
  if (!is.null(by)) nd$.grp <- do.call(paste, c(nd[by], sep = " / "))

  tibble::as_tibble(nd)
}

#' Population-level posterior survival probabilities at a prediction grid
#'
#' Thin wrapper around [brms::posterior_linpred()] with `transform = TRUE` so
#' it returns survival probabilities on the natural scale. Random effects are
#' marginalised out by default (`re_formula = NA`) so the result is a
#' population-level prediction.
#'
#' @param workflow   A fitted `bayes_tls`.
#' @param newdata    Prediction grid from [new_tdt_grid()].
#' @param ndraws     Integer number of posterior draws to use, or `NULL` for all.
#' @param re_formula Passed to [brms::posterior_linpred()]. `NA` (default)
#'                   marginalises random effects.
#' @return Numeric matrix `[ndraws x nrow(newdata)]` of survival probabilities.
#' @keywords internal
posterior_linpred_tdt <- function(workflow,
                                  newdata,
                                  ndraws     = NULL,
                                  re_formula = NA) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  args <- list(object = workflow$fit, newdata = newdata,
               re_formula = re_formula, transform = TRUE)
  if (!is.null(ndraws)) args$ndraws <- ndraws

  do.call(brms::posterior_linpred, args)
}

#' Summarise observed survival counts to mean ± SE per (temp, duration) cell
#'
#' Useful as an overlay on plotted posterior curves. Returns one row per
#' temperature × duration combination, with the mean observed survival
#' proportion across replicates and a standard error.
#'
#' @param observed A tibble from [standardize_data()].
#' @return A tibble with columns `temp`, `duration`, `survival_mean`,
#'         `survival_se`, `survival_lower`, `survival_upper`, `n_units`,
#'         `n_total_sum`.
#' @examples
#' raw <- data.frame(T = rep(30, 6), hrs = rep(c(1, 5), each = 3),
#'                   n = 30, alive = c(28, 27, 29, 10, 12, 8))
#' d <- standardize_data(raw, "T", "hrs", n_total = "n", n_surv = "alive")
#' summarise_observed_survival(d)
#' @export
summarise_observed_survival <- function(observed) {
  # n_total_sum is only available for count data; a continuous-proportion
  # (Beta) frame has no n_total column, so omit that summary there.
  has_n_total <- "n_total" %in% names(observed)
  out <- observed |>
    dplyr::group_by(temp, duration) |>
    dplyr::summarise(
      survival_mean = mean(survival, na.rm = TRUE),
      survival_se   = stats::sd(survival, na.rm = TRUE) / sqrt(dplyr::n()),
      n_units       = dplyr::n(),
      n_total_sum   = if (has_n_total) sum(n_total, na.rm = TRUE)
                      else NA_integer_,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      survival_se    = dplyr::if_else(is.finite(survival_se), survival_se, 0),
      survival_lower = pmax(0, survival_mean - survival_se),
      survival_upper = pmin(1, survival_mean + survival_se)
    )
  if (!has_n_total) out$n_total_sum <- NULL
  out
}

#' Posterior survival curves on a temperature × duration grid
#'
#' Predicts survival probability at every combination of `temps` × `durations`
#' under a fitted 4PL, then returns the posterior median and 95% credible
#' interval at each grid point. Random effects are marginalised out.
#'
#' Use this for *curves* (a handful of temperatures, dense in duration). For a
#' 2-D survival heatmap call [derive_tdt_landscape()] instead.
#'
#' @param workflow  A fitted `bayes_tls`.
#' @param temps     Numeric vector of temperatures (°C). Default: unique assay
#'                  temperatures in the training data.
#' @param durations Numeric vector of durations. Default: 250 log-spaced values
#'                  spanning 0.2× to 5× the observed range.
#' @param ndraws    Integer number of posterior draws to use. Default 1000.
#' @param probs     Numeric length-3 quantile probabilities (lower, median,
#'                  upper). Default `c(0.025, 0.5, 0.975)`.
#' @param by        Optional moderator column(s) for per-group curves. `NULL`
#'                  (default) uses the fit's moderators; a single-condition fit
#'                  returns the ungrouped curves, a grouped fit one block per
#'                  group with the moderator column(s) prepended to `summary`.
#' @return A list with elements `summary` (tibble of `temp`, `duration`,
#'         `survival_median`, `survival_lower`, `survival_upper`; plus the
#'         moderator column(s) for a grouped fit) and `draws`
#'         (the raw posterior matrix as a long tibble).
#' @examples
#' \dontrun{
#' wf <- fit_4pl(d, ...)
#' pred <- predict_survival_curves(wf, temps = c(30, 32, 34))
#' }
#' @export
predict_survival_curves <- function(workflow,
                                    temps     = NULL,
                                    durations = NULL,
                                    ndraws    = 1000,
                                    probs     = c(0.025, 0.5, 0.975),
                                    by        = NULL) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  # Grouped fits predict per group: the moderator column(s) are crossed into the
  # grid so posterior_linpred() validates. `by = NULL` uses the fit's moderators
  # (a single-condition fit has none -> the original ungrouped curves).
  by   <- tdt_resolve_by(workflow, by)
  data <- workflow$data
  if (is.null(temps))     temps     <- sort(unique(data$temp))
  if (is.null(durations)) {
    drange <- range(data$duration, na.rm = TRUE)
    durations <- 10 ^ seq(log10(drange[1] / 5), log10(drange[2] * 5),
                          length.out = 250)
  }

  nd   <- new_tdt_grid(workflow, temps = temps, durations = durations, by = by)
  pred <- posterior_linpred_tdt(workflow, nd, ndraws = ndraws, re_formula = NA)

  summary <- tibble::tibble(
    temp            = nd$temp,
    duration        = nd$duration,
    survival_lower  = apply(pred, 2, stats::quantile, probs = probs[1], na.rm = TRUE),
    survival_median = apply(pred, 2, stats::quantile, probs = probs[2], na.rm = TRUE),
    survival_upper  = apply(pred, 2, stats::quantile, probs = probs[3], na.rm = TRUE)
  )
  if (!is.null(by))   # prepend the moderator column(s) for a grouped fit
    summary <- tibble::as_tibble(cbind(nd[, by, drop = FALSE], summary))

  list(summary = summary, draws_matrix = pred, grid = nd)
}

#' Per-draw threshold inversion of a posterior-prediction matrix
#'
#' Internal helper. For each row (draw) of `pred_mat`, finds the value of `x`
#' at which the predicted probability crosses `target`, via linear
#' interpolation. Returns one value per draw.
#'
#' `stats::approx()` requires its `x` argument (here, the survival values) to
#' be monotonically increasing. Because the 4PL is monotonically *decreasing*
#' in duration (and in temperature for a fixed duration), we sort by `p`
#' ascending before calling `approx()`.
#'
#' @param pred_mat Numeric matrix `[ndraws × ngrid]` of survival probabilities.
#' @param x        Numeric vector of grid values along which to invert
#'                 (length = `ncol(pred_mat)`).
#' @param target   Numeric scalar: the probability value to invert at.
#' @return Numeric vector of length `nrow(pred_mat)`.
#' @keywords internal
threshold_x_by_draw <- function(pred_mat, x, target) {
  apply(pred_mat, 1, function(p) {
    p <- as.numeric(p)
    if (any(!is.finite(p))) return(NA_real_)
    if (target > max(p) || target < min(p)) return(NA_real_)
    ord <- order(p)
    # Duplicate p values at the 4PL asymptote plateaus are mathematically
    # harmless (y is constant there), but approx() warns on every call.
    suppressWarnings(stats::approx(p[ord], x[ord], xout = target)$y)
  })
}
