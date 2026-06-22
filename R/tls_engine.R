# Shared TDT extraction engine. Every post-fit quantity (z, CTmax, T_crit, LT
# curves, survival, heat injury) is derived from brms::posterior_linpred(nlpar=)
# evaluated at an explicit grid — NEVER from raw b_*_Intercept coefficient names.
# posterior_linpred(nlpar = "mid") returns the linear predictor for a midpoint
# fit and the nlf-derived mid for a direct fit, so this engine is identically
# parameterisation- and coding-agnostic. tls() and the single-condition readers
# (extract_tdt, derive_z, derive_temperature_for_duration, tdt_parameter_table,
# extract_4pl_pars) all build on these helpers; the pure-math layer
# (tls_local_z / tls_invert_logLT) operates on already-evaluated matrices so the
# z/CTmax math stays fast-testable without a Stan fit.

# Build a moderator x temperature grid for posterior_linpred. Crosses `by`-levels
# with `temp_grid`, fills any model-data columns the grid lacks (so brms
# validates), and tags groups in `.grp` (`"all"` when `by = NULL`). `temp` is the
# (centred) temperature column name in the model data.
tls_build_grid <- function(mdata, by = NULL, temp = "temp_c", temp_grid = NULL,
                           newdata = NULL) {
  if (is.null(newdata)) {
    if (is.null(temp_grid))
      temp_grid <- seq(min(mdata[[temp]]), max(mdata[[temp]]), length.out = 11)
    base <- if (is.null(by)) data.frame(.dummy = 1)
            else unique(mdata[, by, drop = FALSE])
    newdata <- do.call(rbind, lapply(seq_len(nrow(base)), function(i) {
      g <- base[rep(i, length(temp_grid)), , drop = FALSE]
      g[[temp]] <- temp_grid
      g
    }))
    newdata$.dummy <- NULL
    rownames(newdata) <- NULL
  }
  for (v in setdiff(names(mdata), names(newdata))) newdata[[v]] <- mdata[[v]][1]
  newdata$.grp <- if (is.null(by)) "all" else
    do.call(paste, c(newdata[by], sep = " / "))
  newdata
}

# Resolve the grouping columns for a group-aware reader: an explicit `by` wins;
# otherwise use the fit's moderators (`meta$group_vars`); a single-condition fit
# has none -> NULL -> one group "all" -> the reader's single-condition output.
# Safety net: a fit that LOOKS grouped (multi-level CTmax/z or mid coefficients)
# but records no moderators (e.g. a hand-wrapped brmsfit) must not be silently
# extracted at its reference level (tls_build_grid would fill the moderator from
# row 1) — require the caller to name the moderator(s).
tdt_resolve_by <- function(workflow, by) {
  if (!is.null(by)) return(by)
  gv <- workflow$meta$group_vars
  if (length(gv)) return(gv)
  if (tdt_is_grouped(workflow))
    stop("This fit appears to vary CTmax/z (or mid) by a moderator, but its ",
         "moderator column(s) are not recorded in meta$group_vars. Pass ",
         "`by = <column(s)>` to extract per group (or use tls(object, by = ...)).",
         call. = FALSE)
  NULL
}

# Subsample draw indices reproducibly, matching the retired tdt_extract_pars()
# selection (`d[sort(sample.int(n, ndraws)), ]`) so the single-condition readers
# stay bit-identical to the coefficient-parsing implementation. Returns NULL
# (use all draws, in order) when `ndraws` is NULL or >= the available draws. The
# caller is responsible for set.seed() beforehand.
tls_draw_ids <- function(fit, ndraws) {
  if (is.null(ndraws)) return(NULL)
  n <- brms::ndraws(fit)
  if (!is.finite(ndraws) || ndraws >= n) return(NULL)
  sort(sample.int(n, ndraws))
}

# Evaluate the four 4PL sub-parameters at `newdata`, returning per-draw matrices
# (`[ndraws x nrow(newdata)]`) on the natural scale plus the threshold-specific
# logLT (relative: mid; absolute p: mid + log((up - p)/(p - low))/k). This is the
# ONLY engine layer that calls brms — the gate boundary for fast vs gated tests.
# `draw_ids` selects specific posterior draws (see tls_draw_ids); `ndraws` is the
# brms count-based subsample (used by tls()); pass at most one.
tls_eval_subpars <- function(fit, newdata, bounds,
                             nlpars = c("lowraw", "upraw", "logk", "mid"),
                             re_formula = NA, ndraws = NULL, draw_ids = NULL,
                             mode = "absolute", p = 0.5) {
  lp <- function(np) brms::posterior_linpred(fit, newdata = newdata, nlpar = np,
                                             re_formula = re_formula, ndraws = ndraws,
                                             draw_ids = draw_ids)
  low <- bounds$low_min + stats::plogis(lp(nlpars[1])) * bounds$low_w
  up  <- bounds$up_min  + stats::plogis(lp(nlpars[2])) * bounds$up_w
  k   <- exp(lp(nlpars[3]))
  mid <- lp(nlpars[4])
  logLT <- if (identical(mode, "absolute")) mid + log((up - p) / (p - low)) / k else mid
  list(low = low, up = up, k = k, mid = mid, logLT = logLT)
}

# --- pure math on already-evaluated matrices (no brms; fast-testable) ---------

# Local z(T) = -1 / d(logLT)/dT by central finite difference, plus the pooled z
# (per-draw mean over the grid). `logLT_plus`/`logLT_minus` are the logLT
# matrices evaluated at `temp_grid + h` / `temp_grid - h`. For a LINEAR logLT
# (relative threshold, midpoint or direct) the central difference is EXACT, so
# pooled z is bit-identical to the closed-form -1/b_mid_temp_c. Mirrors the math
# of the retired tdt_z_from_pars(), consuming matrices instead of coefficients.
tls_local_z <- function(logLT_plus, logLT_minus, h, temp_grid,
                        probs = c(0.025, 0.5, 0.975), local = TRUE) {
  np    <- nrow(logLT_plus)
  slope <- (logLT_plus - logLT_minus) / (2 * h)
  zloc  <- -1 / slope
  zloc[!is.finite(zloc)] <- NA_real_
  z_pooled <- rowMeans(zloc, na.rm = TRUE)
  z_pooled[is.nan(z_pooled)] <- NA_real_

  draws <- tibble::tibble(.draw = seq_len(np), z = z_pooled) |>
    dplyr::filter(is.finite(z))
  qz <- stats::quantile(draws$z, probs, names = FALSE, na.rm = TRUE)
  if (isTRUE(local)) {
    local_draws <- tibble::tibble(
      .draw = rep(seq_len(np), times = length(temp_grid)),
      temp  = rep(temp_grid, each = np),
      z     = as.vector(zloc)
    ) |> dplyr::filter(is.finite(z))
    local_summary <- local_draws |>
      dplyr::group_by(temp) |>
      dplyr::summarise(z_median = stats::median(z),
                       z_lower  = stats::quantile(z, probs[1], names = FALSE),
                       z_upper  = stats::quantile(z, probs[3], names = FALSE),
                       .groups  = "drop")
  } else {
    local_draws <- NULL; local_summary <- NULL
  }
  list(draws = draws,
       summary = tibble::tibble(z_median = qz[2], z_lower = qz[1], z_upper = qz[3]),
       local_draws = local_draws, local_summary = local_summary)
}

# Per-draw CTmax: temperature where logLT(T) = target, by vectorised inverse
# linear interpolation across draws (the curve is monotone-decreasing in T for
# the disjoint-bounds 4PL), with an exact per-row stats::approx fallback for any
# non-monotone / non-finite draw. `M` is the evaluated logLT matrix
# `[ndraws x length(temp_grid)]`. Mirrors the retired tdt_ctmax_from_pars().
tls_invert_logLT <- function(M, target, temp_grid) {
  np <- nrow(M); nT <- length(temp_grid)
  Tc <- rep(NA_real_, np)
  j  <- rowSums(M > target, na.rm = TRUE)
  ok <- is.finite(j) & j >= 1L & j <= (nT - 1L)
  iv <- which(ok); jj <- j[iv]
  y1 <- M[cbind(iv, jj)]; y2 <- M[cbind(iv, jj + 1L)]
  Tc[iv] <- temp_grid[jj] +
    (target - y1) * (temp_grid[jj + 1L] - temp_grid[jj]) / (y2 - y1)
  dM  <- M[, -1, drop = FALSE] - M[, -nT, drop = FALSE]
  bad <- !is.finite(rowSums(M)) | rowSums(dM >= 0, na.rm = TRUE) > 0
  for (i in which(bad)) {
    y <- M[i, ]; fin <- is.finite(y)
    if (sum(fin) < 2L) { Tc[i] <- NA_real_; next }
    o <- order(y[fin])
    Tc[i] <- suppressWarnings(stats::approx(y[fin][o], temp_grid[fin][o],
                                            xout = target)$y)
  }
  Tc
}
