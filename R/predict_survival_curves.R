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
#' @param workflow  A `tdt_4pl_workflow` object.
#' @param temps     Numeric vector of assay temperatures (°C).
#' @param durations Numeric vector of exposure durations (same unit as the
#'                  training data, typically hours).
#' @param n_total   Integer trials count. `NULL` (default) uses the median of
#'                  the training data.
#' @return A tibble with one row per `(temp, duration)` pair.
#' @keywords internal
new_tdt_grid <- function(workflow,
                         temps,
                         durations,
                         n_total = NULL) {
  data <- workflow$data
  meta <- workflow$meta

  if (is.null(n_total)) {
    n_total <- as.integer(round(stats::median(data$n_total, na.rm = TRUE)))
  }

  nd <- expand.grid(temp     = temps,
                    duration = durations,
                    KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)
  nd$n_total <- n_total
  nd$logd    <- log10(nd$duration)
  nd$temp_c  <- nd$temp - meta$temp_mean

  for (re_var in tdt_random_effect_variables(meta$random_effects)) {
    if (re_var %in% names(data)) {
      if (is.factor(data[[re_var]])) {
        nd[[re_var]] <- factor(levels(data[[re_var]])[1],
                               levels = levels(data[[re_var]]))
      } else {
        nd[[re_var]] <- data[[re_var]][1]
      }
    }
  }

  tibble::as_tibble(nd)
}

#' Population-level posterior survival probabilities at a prediction grid
#'
#' Thin wrapper around [brms::posterior_linpred()] with `transform = TRUE` so
#' it returns survival probabilities on the natural scale. Random effects are
#' marginalised out by default (`re_formula = NA`) so the result is a
#' population-level prediction.
#'
#' @param workflow   A fitted `tdt_4pl_workflow`.
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
  observed |>
    dplyr::group_by(temp, duration) |>
    dplyr::summarise(
      survival_mean = mean(survival, na.rm = TRUE),
      survival_se   = stats::sd(survival, na.rm = TRUE) / sqrt(dplyr::n()),
      n_units       = dplyr::n(),
      n_total_sum   = sum(n_total, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      survival_se    = dplyr::if_else(is.finite(survival_se), survival_se, 0),
      survival_lower = pmax(0, survival_mean - survival_se),
      survival_upper = pmin(1, survival_mean + survival_se)
    )
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
#' @param workflow  A fitted `tdt_4pl_workflow`.
#' @param temps     Numeric vector of temperatures (°C). Default: unique assay
#'                  temperatures in the training data.
#' @param durations Numeric vector of durations. Default: 250 log-spaced values
#'                  spanning 0.2× to 5× the observed range.
#' @param ndraws    Integer number of posterior draws to use. Default 1000.
#' @param probs     Numeric length-3 quantile probabilities (lower, median,
#'                  upper). Default `c(0.025, 0.5, 0.975)`.
#' @return A list with elements `summary` (tibble of `temp`, `duration`,
#'         `survival_median`, `survival_lower`, `survival_upper`) and `draws`
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
                                    probs     = c(0.025, 0.5, 0.975)) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  data <- workflow$data
  if (is.null(temps))     temps     <- sort(unique(data$temp))
  if (is.null(durations)) {
    drange <- range(data$duration, na.rm = TRUE)
    durations <- 10 ^ seq(log10(drange[1] / 5), log10(drange[2] * 5),
                          length.out = 250)
  }

  nd   <- new_tdt_grid(workflow, temps = temps, durations = durations)
  pred <- posterior_linpred_tdt(workflow, nd, ndraws = ndraws, re_formula = NA)

  summary <- tibble::tibble(
    temp            = nd$temp,
    duration        = nd$duration,
    survival_lower  = apply(pred, 2, stats::quantile, probs = probs[1], na.rm = TRUE),
    survival_median = apply(pred, 2, stats::quantile, probs = probs[2], na.rm = TRUE),
    survival_upper  = apply(pred, 2, stats::quantile, probs = probs[3], na.rm = TRUE)
  )

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
    stats::approx(p[ord], x[ord], xout = target)$y
  })
}
