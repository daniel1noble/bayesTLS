# Heat injury and predicted survival under a fluctuating temperature trace.
# Carries the model's full posterior through to HI(t) and S(t) trajectories
# with credible intervals. The algorithm:
#
#   1. Extract per-draw 4PL parameters (low, up, k, mid_int, mid_temp).
#   2. At each time step, compute LT_target(T) analytically from the 4PL,
#      giving damage rate = 1 / LT_target.
#   3. Optionally add a temperature-dependent repair rate (Sharpe-Schoolfield).
#   4. Forward-Euler integrate (damage - repair) across the trace, using each
#      interval's own width (dose starts at zero at the first time point).
#   5. Map cumulative dose back to predicted survival via the fitted model's
#      curve (log-logistic, log-normal or Weibull; see R/models.R). With
#      shape = "varying", survival follows each individual's injury against
#      its own threshold time at every temperature, or the cumulative-exposure
#      update (see hi_rank_integrate()).

#' Per-draw natural-scale 4PL parameters from a fitted workflow
#'
#' Faithfully extracts ALL FOUR natural-scale 4PL sub-parameters
#' (`low`, `up`, `k`, `mid`) per posterior draw, at one or more assay
#' temperatures and per moderator group, retaining whatever temperature and/or
#' group structure each sub-parameter was fit with (`~ 1`, `~ temp_c`,
#' `~ temp_c * group`, random effects, …). It surfaces the shared TDT engine
#' (`tls_eval_subpars()` over a `tls_build_grid()` temperature × moderator
#' grid via `brms::posterior_linpred(nlpar=)`), so the result is
#' parameterisation- and coding-agnostic — no coefficient-name parsing.
#'
#' `mid` is the natural-scale midpoint sub-parameter (the log10 time at the
#' per-draw relative threshold) evaluated at each temperature; on an absolute
#' fit it is still the bare `mid` nlpar (the asymmetry correction is applied
#' downstream by the threshold-specific readers, not folded in here).
#'
#' This replaces the earlier constant-shape extractor that discarded the
#' temperature slopes on `low`/`up`/`k`. The constant-in-T view the
#' heat-injury integral needs (`low, up, k` at the centring temperature plus
#' a linear `mid_int`/`mid_temp`) now lives in the internal helper
#' [hi_pars()], which [predict_heat_injury()] calls for `shape = "constant"`.
#'
#' @param workflow Fitted `bayes_tls`.
#' @param temps Numeric vector of assay temperatures (°C, uncentred) to
#'   evaluate the curve parameters at. `NULL` (default) uses the fit's observed
#'   unique assay temperatures.
#' @param by Optional moderator column(s) for per-group parameters. `NULL`
#'   (default) uses the fit's moderators (all levels); a single-condition fit
#'   returns the ungrouped tibble. A grouped fit prepends the moderator
#'   column(s).
#' @param re_formula Passed to [brms::posterior_linpred()]. `NA` (default)
#'   marginalises out group-level (random) effects so the result is a
#'   population-level curve; `NULL` conditions on the fitted random effects.
#' @param ndraws Posterior draws to use; `NULL` (default) uses the full
#'   posterior. Pass an integer to subsample (the caller is responsible for
#'   `set.seed()` if a reproducible subsample is wanted).
#' @return A long tibble with `(.draw, temp, low, up, k, mid)` columns (plus
#'         the moderator column(s) for a grouped fit, immediately after
#'         `.draw`), one row per draw × temperature × group, filtered to rows
#'         producing valid parameter values (`k > 0`, `up > low`, all finite).
#'         `temp` is the uncentred assay temperature (°C).
#' @seealso [hi_pars()] for the constant-shape view used by the heat-injury
#'   integral; [get_4pl_est()] which wraps this for the `"draws"` view.
#' @examples
#' \dontrun{
#' pars <- extract_4pl_pars(wf)                 # observed temps, fit's groups
#' extract_4pl_pars(wf, temps = c(30, 33, 36))  # at chosen temperatures
#' head(pars)
#' }
#' @export
extract_4pl_pars <- function(workflow, temps = NULL, by = NULL,
                             re_formula = NA, ndraws = NULL) {
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)

  fit       <- get_brmsfit(workflow)
  by        <- tdt_resolve_by(workflow, by)
  temp_mean <- workflow$meta$temp_mean %||% 0

  # Default to the fit's observed unique assay temperatures (uncentred °C).
  if (is.null(temps)) {
    temps <- sort(unique(fit$data$temp_c)) + temp_mean
  }
  temp_grid <- temps - temp_mean                       # centre for the engine

  draw_ids <- if (is.null(ndraws)) NULL else tls_draw_ids(fit, ndraws)
  nd <- tls_build_grid(fit$data, by = by, temp = "temp_c", temp_grid = temp_grid)
  sp <- tls_eval_subpars(fit, nd, tdt_workflow_bounds(workflow), mode = "relative",
                         re_formula = re_formula, draw_ids = draw_ids)

  ndr <- nrow(sp$mid)
  per_group <- lapply(unique(nd$.grp), function(g) {
    gi <- which(nd$.grp == g)
    # One block per (draw, temperature) cell in this group, in temp order.
    blocks <- lapply(gi, function(col) {
      tibble::tibble(
        .draw = seq_len(ndr),
        temp  = nd$temp_c[col] + temp_mean,
        low   = sp$low[, col], up = sp$up[, col], k = sp$k[, col],
        mid   = sp$mid[, col]
      )
    })
    out <- dplyr::bind_rows(blocks)
    out <- dplyr::filter(out,
                         is.finite(low) & is.finite(up) & is.finite(k) &
                         is.finite(mid) & k > 0 & up > low)
    if (!is.null(by)) {
      gcols <- nd[gi[1], by, drop = FALSE]
      out <- cbind(gcols[rep(1, nrow(out)), , drop = FALSE], out, row.names = NULL)
    }
    out
  })
  tibble::as_tibble(dplyr::bind_rows(per_group))
}

#' Constant-shape 4PL parameters for the heat-injury integral (internal)
#'
#' The classical heat-injury integral is evaluated under the assumption that
#' the asymptotes (`low`, `up`) and steepness (`k`) are constant in
#' temperature and only the midpoint shifts. This helper returns exactly that
#' view, per posterior draw and per moderator group: `low`, `up`, `k` read at
#' the centring temperature (`temp_c = 0`, i.e. `T_bar`), plus a linear
#' midpoint decomposed into intercept (`mid_int`, at `T_bar`) and slope
#' (`mid_temp`, per °C). This is the per-draw input
#' [predict_heat_injury()] uses for `shape = "constant"`.
#'
#' Reads `posterior_linpred(nlpar=)` at `temp_c = 0` (asymptotes + midpoint
#' intercept) and `temp_c = 1` (for the linear midpoint slope), so it is
#' parameterisation- and coding-agnostic. If the temperature slopes on
#' `low`/`up`/`k` are shrunk near zero by the data this introduces no bias;
#' otherwise it is the same approximation the classical HI framework makes
#' (see `shape = "varying"` for the slope-aware alternative).
#'
#' @param workflow Fitted `bayes_tls`.
#' @param by Resolved moderator column(s) (already through `tdt_resolve_by()`),
#'   or `NULL` for a single-condition fit.
#' @return A tibble with `(.draw, low, up, k, mid_int, mid_temp)` (plus the
#'   moderator column(s) for a grouped fit), filtered to valid draws.
#' @keywords internal
hi_pars <- function(workflow, by = NULL) {
  fit <- get_brmsfit(workflow)
  nd  <- tls_build_grid(fit$data, by = by, temp = "temp_c", temp_grid = c(0, 1))
  sp  <- tls_eval_subpars(fit, nd, tdt_workflow_bounds(workflow), mode = "relative")

  per_group <- lapply(unique(nd$.grp), function(g) {
    gi <- which(nd$.grp == g)
    c0 <- gi[which(nd$temp_c[gi] == 0)]; c1 <- gi[which(nd$temp_c[gi] == 1)]
    out <- tibble::tibble(
      .draw    = seq_len(nrow(sp$mid)),
      low      = sp$low[, c0], up = sp$up[, c0], k = sp$k[, c0],
      mid_int  = sp$mid[, c0], mid_temp = sp$mid[, c1] - sp$mid[, c0]
    )
    out <- dplyr::filter(out,
                         is.finite(low) & is.finite(up) & is.finite(k) &
                         is.finite(mid_int) & is.finite(mid_temp) & k > 0 & up > low)
    if (!is.null(by)) out <- cbind(nd[gi[1], by, drop = FALSE], out, row.names = NULL)
    out
  })
  dplyr::bind_rows(per_group)
}

#' Temperature-local 4PL curves per posterior draw (internal, for shape="varying")
#'
#' Evaluates the full natural-scale 4PL `low(T)`, `up(T)`, `k(T)`, `mid(T)` per
#' posterior draw at a set of unique temperatures, retaining every sub-parameter's
#' temperature/group structure (the slopes the constant-shape [hi_pars()] view
#' drops). Returns per-draw × per-temperature matrices, keyed for fast lookup by
#' the shape-varying heat-injury integrator: each trace point's temperature maps
#' to a column, each posterior draw to a row.
#'
#' @param workflow Fitted `bayes_tls`.
#' @param by Resolved moderator column(s), or `NULL` for a single-condition fit.
#' @param temps_unique Numeric vector of unique (uncentred °C) temperatures to
#'   evaluate the curves at -- typically `unique(trace$temp)`.
#' @return A named list, one element per moderator group (named by the group
#'   label, `"all"` when ungrouped). Each element is a list with:
#'   `temps` (the input temperatures, in evaluation order), `.draw` (the
#'   posterior draw indices kept after filtering invalid draws), and the
#'   matrices `low`, `up`, `k`, `mid`, each `[length(.draw) × length(temps)]`.
#'   A draw is kept only if ALL its cells are valid (finite, `k > 0`,
#'   `up > low`) at every temperature, so the integrator never hits a partial
#'   draw.
#' @keywords internal
hi_local_curves <- function(workflow, by = NULL, temps_unique) {
  fit       <- get_brmsfit(workflow)
  temp_mean <- workflow$meta$temp_mean %||% 0
  temp_grid <- temps_unique - temp_mean
  nd <- tls_build_grid(fit$data, by = by, temp = "temp_c", temp_grid = temp_grid)
  sp <- tls_eval_subpars(fit, nd, tdt_workflow_bounds(workflow), mode = "relative")
  ndr <- nrow(sp$mid)

  out <- lapply(unique(nd$.grp), function(g) {
    gi <- which(nd$.grp == g)
    # Columns of nd for this group, ordered to match temps_unique.
    ord <- gi[order(match(round(nd$temp_c[gi] + temp_mean, 6),
                          round(temps_unique, 6)))]
    low <- sp$low[, ord, drop = FALSE]; up <- sp$up[, ord, drop = FALSE]
    k   <- sp$k[, ord, drop = FALSE];  mid <- sp$mid[, ord, drop = FALSE]
    # Keep only draws valid at EVERY temperature (so no partial-draw traces).
    ok <- apply(is.finite(low) & is.finite(up) & is.finite(k) & is.finite(mid) &
                  (k > 0) & (up > low), 1, all)
    list(temps = temps_unique, .draw = which(ok), grp = g,
         low = low[ok, , drop = FALSE], up = up[ok, , drop = FALSE],
         k = k[ok, , drop = FALSE],   mid = mid[ok, , drop = FALSE],
         by_row = if (!is.null(by)) nd[gi[1], by, drop = FALSE] else NULL)
  })
  names(out) <- vapply(out, function(x) x$grp, character(1))
  out
}

#' Analytical inverse 4PL: duration to reach a target survival at a given temperature
#'
#' For each `temp`, computes the exposure duration at which the 4PL gives
#' the chosen survival threshold. In `"relative"` mode the threshold is the
#' per-draw midpoint between the asymptotes, which collapses to the bare 4PL
#' `mid` parameter on the natural time axis. In absolute mode it is the
#' literal probability `survival_target`.
#'
#' @param temp           Numeric vector of temperatures (°C).
#' @param survival_target Either the literal probability to invert at (must
#'                       lie strictly between `low` and `up`), or the string
#'                       `"relative"` to use `(low + up)/2`.
#' @param low,up,k       Scalar (or per-temp) 4PL parameters.
#' @param mid_int,mid_temp Midpoint sub-model coefficients.
#' @param temp_mean      Centring temperature used by the model.
#' @param model          The dose-response model: a name (default
#'                       `"log-logistic"`) or a resolved registry entry. Only
#'                       the absolute threshold depends on it.
#' @return Numeric vector of durations in the model's time units. `NA` where
#'         the threshold is outside `(low, up)` in absolute mode.
#' @keywords internal
time_to_surv_threshold_4pl <- function(temp, survival_target,
                                       low, up, k,
                                       mid_int, mid_temp, temp_mean,
                                       model = "log-logistic") {
  temp_c <- temp - temp_mean
  mid    <- mid_int + mid_temp * temp_c
  if (is.character(survival_target) && length(survival_target) == 1L &&
      survival_target == "relative") {
    return(10 ^ mid)
  }
  if (survival_target <= low || survival_target >= up)
    return(rep(NA_real_, length(temp)))
  # The model's offset; for the log-logistic log((up - p)/(p - low))/k.
  10 ^ (mid + tdt_as_model_spec(model)$abs_corr(low, up, k, survival_target))
}

#' Survival corresponding to an accumulated dose
#'
#' Inverts the dose interpretation: at `dose = 1`, survival equals the
#' threshold by construction. The 4PL maps dose smoothly to survival on the
#' reference dose-response curve.
#'
#' @param dose Numeric vector of cumulative doses (in units where 1 dose =
#'             1 LT-dose at the chosen threshold).
#' @param low,up,k 4PL parameters at the reference (centring) temperature.
#' @param target_surv Either the literal probability defining "1 dose" or the
#'                    string `"relative"` for the `(low + up)/2` threshold.
#' @param model       The dose-response model: a name (default
#'                    `"log-logistic"`) or a resolved registry entry. Survival
#'                    is read off that model's curve.
#' @return Numeric vector of predicted survival probabilities. `NA` for every
#'         element when a numeric `target_surv` lies outside `(low, up)` and is
#'         therefore unreachable on this draw's 4PL (mirrors
#'         [time_to_surv_threshold_4pl()], so the draw drops cleanly rather than
#'         producing `log()`-of-negative `NaN`s that would poison the trace).
#' @keywords internal
survival_from_dose <- function(dose, low, up, k, target_surv = "relative",
                               model = "log-logistic") {
  mod      <- tdt_as_model_spec(model)
  dose_use <- pmax(dose, 1e-12)
  if (is.character(target_surv) && length(target_surv) == 1L &&
      target_surv == "relative") {
    c_target <- 0
  } else {
    # Unreachable target on this draw -- guard like time_to_surv_threshold_4pl()
    # so the whole draw drops as NA instead of evaluating the correction outside
    # its domain (which would yield NaN and silently bias the medians upward).
    if (target_surv <= low || target_surv >= up)
      return(rep(NA_real_, length(dose)))
    c_target <- mod$abs_corr(low, up, k, target_surv)
  }
  # At dose 1 the curve sits exactly at the threshold: u = k * c_target.
  mod$curve(low, up, k * (log10(dose_use) + c_target))
}

#' Varying-shape survival from individual thresholds (internal)
#'
#' The integrator behind `predict_heat_injury(shape = "varying")`. An
#' individual of survival rank `p` is the one that dies, at a constant
#' temperature `T`, when population survival there falls to `p`: at the
#' threshold time `LT_p(T)`. Each individual accumulates injury additively,
#' `sum(dt / LT_p(T))` (the TDT additivity rule applied to the individual), and
#' dies when it reaches 1. Injury is kept as relative dose in bins by the
#' temperature at which it was accrued, so individual `p`'s injury is
#' `sum_T rho_p(T) A_T` with `rho_p(T) = LT_rel(T) / LT_p(T) = 10^(-c_p(T))`,
#' `c_p` being the model's absolute offset (`abs_corr`). Repair removes the
#' same fraction from every bin, so the total relative dose falls by the repair
#' amount, as in the constant shape. Because `rho_p(T)` increases with `p`,
#' the survivors are the ranks below the one whose injury equals 1; survival
#' is that rank, as a running minimum when mortality is irreversible.
#' At a constant temperature this is the fitted curve, and with `low`, `up`
#' and `k` the same at every temperature (and `mid` linear in it) it is the
#' constant shape exactly.
#'
#' Survival starts at the top of the curve at the first temperature where the
#' trace accrues injury (above `T_c`, over a step of positive width); a trace
#' that accrues none keeps `s_none`. Ranks at or above `up(T)` die on exposure
#' to `T`; those deaths are permanent even when mortality is reversible, as the
#' share `1 - up` is in the constant shape.
#'
#' @param L,U,K,MID Matrices `[bins x draws]` of `low`, `up`, `k` and `mid`
#'   (log10 model time units) at each bin's temperature.
#' @param bin Integer bin of each trace point.
#' @param dt Interval widths (hours), length `length(bin) - 1`.
#' @param below Logical per trace point: at or below `T_c` (no injury).
#' @param unit_h Hours per model time unit.
#' @param mod Resolved model registry entry.
#' @param p_target `NULL` for heat injury in relative doses, or an absolute
#'   survival target whose rank's injury is reported.
#' @param rep_rate `NULL`, or the repair rate at each trace point (doses of
#'   the target per hour).
#' @param conv Per-draw factor converting `rep_rate` to relative doses.
#' @param scale_with_surv Scale repair by survival / `up` at the starting
#'   temperature (at most 1).
#' @param irreversible Survival cannot rebound.
#' @param s_none Survival, per draw, when no step accrues injury (the constant
#'   shape's starting value); `NULL` uses the first trace temperature's top.
#' @param method `"individual"` (above) or `"cumulative"`: cumulative exposure,
#'   where survival is the state. A step that injures moves each draw's state
#'   to the step's curve, at the dose where that curve has the state's
#'   survival, and adds the step's dose; a state at or below that curve's floor
#'   stays on its own curve, and the step adds nothing to it. At or below `T_c`
#'   the state stays where it is. Repair lowers the dose on the curve the state
#'   lies on. With irreversible mortality survival is the lowest survival of
#'   the injury so far. Heat injury is the same for both unless repair scales
#'   with survival.
#' @return A list with `survival` and `dose`, matrices `[points x draws]`, and
#'   `max_iter`, the most solver steps any root of the per-individual method
#'   took.
#' @keywords internal
hi_rank_integrate <- function(L, U, K, MID, bin, dt, below, unit_h, mod,
                              p_target = NULL, rep_rate = NULL, conv = 1,
                              scale_with_surv = TRUE, irreversible = TRUE,
                              s_none = NULL,
                              method = c("individual", "cumulative")) {
  method <- match.arg(method)
  nb <- nrow(L); nd <- ncol(L); n <- length(bin)
  # Steps that accrue injury: of positive width, above T_c. With none, nothing
  # happens: survival keeps its starting value and heat injury stays 0.
  inj <- which(dt > 0 & !below[-n])
  if (!length(inj)) {
    s0 <- if (is.null(s_none)) {
      mod$curve(L[bin[1], ], U[bin[1], ], K[bin[1], ] * log10(1e-12))
    } else s_none
    return(list(survival = matrix(s0, n, nd, byrow = TRUE), dose = matrix(0, n, nd),
                max_iter = 0L))
  }
  ln10 <- log(10); eps <- .Machine$double.eps
  A <- matrix(0, nb, nd)                     # outstanding relative dose by bin
  rows <- integer(0)                         # bins that have received injury
  max_iter <- 0L                             # most Newton/bisection steps in a solve
  lo <- hi <- rep(Inf, nd)                   # lowest low / up over bins holding injury
  # Injury of rank p[i] for draw cols[i] (over the bins holding injury), and
  # optionally its derivative in p.
  rank_injury <- function(p, cols = seq_len(nd), deriv = FALSE) {
    nc <- length(cols)
    if (!length(rows)) {
      return(if (deriv) list(val = numeric(nc), dval = numeric(nc)) else numeric(nc))
    }
    Ac <- A[rows, cols, drop = FALSE]; Lc <- L[rows, cols, drop = FALSE]
    Uc <- U[rows, cols, drop = FALSE]
    P   <- matrix(p, length(rows), nc, byrow = TRUE)
    act <- Ac > 0
    ins <- act & Lc < P & P < Uc
    rho <- matrix(0, length(rows), nc)
    Kc  <- K[rows, cols, drop = FALSE][ins]
    if (any(ins)) rho[ins] <- 10^(-mod$abs_corr(Lc[ins], Uc[ins], Kc, P[ins]))
    val <- colSums(rho * Ac)
    val[colSums(act & P >= Uc) > 0] <- Inf   # the target is already passed
    if (!deriv) return(val)
    drho <- matrix(0, length(rows), nc)
    if (any(ins)) {
      drho[ins] <- -ln10 * rho[ins] * mod$abs_corr_dp(Lc[ins], Uc[ins], Kc, P[ins])
    }
    list(val = val, dval = colSums(drho * Ac))
  }
  # Highest surviving rank, per draw: the p whose injury is exactly 1 (1 when
  # nothing is outstanding). Solved in s, p = hi - (hi - lo) exp(-s), which
  # straightens the singularity at hi (valid curves have low < up, so
  # lo < hi). Newton steps on log injury stay inside a bracket [sa, sb] that
  # holds the root; a step that would leave it, or a Newton step that did not
  # halve the residual, is replaced by bisection. Injury has kinks with an
  # infinite slope where p crosses a bin's `low`, so a short Newton step is
  # not taken as convergence: a root is accepted only when log injury is
  # within 1e-13 of 0 or the bracket has closed, in s or in p (to a few
  # rounding steps). s_max keeps p a few floating-point steps below hi.
  rank_root <- function(guess) {
    root <- rep(1, nd)
    cols <- which(is.finite(lo))
    if (!length(cols)) return(root)
    l <- lo[cols]; h <- hi[cols]; w <- h - l
    s_max <- log(w / (8 * eps * h))
    s  <- pmin(1, s_max / 2)
    g  <- guess[cols]
    ok_g <- is.finite(g) & g > l & g < h     # warm start from the previous root
    s[ok_g] <- pmin(-log((h[ok_g] - g[ok_g]) / w[ok_g]), s_max[ok_g])
    sa <- rep(0, length(cols)); sb <- s_max
    res_prev   <- rep(Inf, length(cols))     # |log injury| where the last step began
    was_newton <- rep(FALSE, length(cols))
    todo <- seq_along(cols)
    it <- 0L
    while (length(todo) && it < 200L) {
      it <- it + 1L
      st <- s[todo]; ht <- h[todo]; wt <- w[todo]; e_s <- exp(-st)
      ev  <- rank_injury(ht - wt * e_s, cols[todo], deriv = TRUE)
      Hs  <- log(ev$val)                     # log injury: the root is at 0
      dHs <- ev$dval / ev$val * wt * e_s
      neg <- Hs < 0
      sa[todo][neg]  <- st[neg]
      sb[todo][!neg] <- st[!neg]
      at <- sa[todo]; bt <- sb[todo]
      done <- abs(Hs) <= 1e-13 | (bt - at) <= 1e-13 * pmax(1, bt) |
        wt * (exp(-at) - exp(-bt)) <= 4 * eps * ht
      s_new <- st - Hs / dHs
      bis <- !is.finite(s_new) | s_new <= at | s_new >= bt |
        (was_newton[todo] & abs(Hs) > 0.5 * res_prev[todo])
      s_new[bis] <- (at[bis] + bt[bis]) / 2
      res_prev[todo]   <- abs(Hs)
      was_newton[todo] <- !bis
      s[todo[!done]] <- s_new[!done]
      todo <- todo[!done]
    }
    max_iter <<- max(max_iter, it)
    root[cols] <- h - w * exp(-s)
    root
  }

  # Cumulative exposure: the state is the curve argument u = k log10(D*) of
  # the current injury, on each draw's curve `cb` (that of the last step that
  # injured it), D* being the relative dose at which that curve has the
  # injury's survival. When a step injures at a new temperature, the state is
  # carried to that step's curve with its survival unchanged, through the log
  # of its deficit from the top or of its excess above the floor, whichever is
  # smaller, so nothing underflows however close survival is to an asymptote.
  # A state at or above the new curve's top starts that curve (deaths on
  # exposure); a state at or below its floor stays on its own curve, and that
  # temperature adds no injury to it.
  log_shift <- function(lx, d) {             # log(exp(lx) + d); -Inf when <= 0
    out <- lx
    pos <- which(d > 0)
    if (length(pos)) {
      ld <- log(d[pos]); x <- lx[pos]; m <- pmax(x, ld)
      out[pos] <- m + log(exp(x - m) + exp(ld - m))
    }
    neg <- which(d < 0)
    if (length(neg)) {
      ld <- log(-d[neg]); x <- lx[neg]; ok <- x > ld
      res <- rep(-Inf, length(neg))
      res[ok] <- x[ok] + log1p(-exp(ld[ok] - x[ok]))
      out[neg] <- res
    }
    out
  }
  ce_map <- function(u, a, b, cols) {        # draws `cols`, from curves `a` to curve b
    ia <- cbind(a, cols)
    La <- L[ia]; Ua <- U[ia]; Lb <- L[b, cols]; Ub <- U[b, cols]
    out  <- u; stay <- rep(FALSE, length(cols))
    todo <- which(!(La == Lb & Ua == Ub))     # same asymptotes: u carries over
    if (!length(todo)) return(list(u = out, stay = stay))
    uu  <- u[todo]
    lra <- log(Ua[todo] - La[todo]); lrb <- log(Ub[todo] - Lb[todo])
    ldf <- log_shift(lra + mod$logGc(uu), Ub[todo] - Ua[todo])   # log(up_b - S)
    lex <- log_shift(lra + mod$logG(uu), La[todo] - Lb[todo])    # log(S - low_b)
    top  <- ldf == -Inf                       # at or above the new top
    bot  <- lex == -Inf & !top                # at or below the new floor: stays
    mid  <- !top & !bot
    by_c <- mid & ldf <= lex
    res  <- uu
    res[top] <- -Inf
    res[by_c] <- mod$Gcinv_log(ldf[by_c] - lrb[by_c])
    res[mid & !by_c] <- mod$Ginv_log(lex[mid & !by_c] - lrb[mid & !by_c])
    out[todo] <- res; stay[todo] <- bot
    list(u = out, stay = stay)
  }
  # Survival of each draw's state on its curve (zero dose read as the constant
  # shape reads it).
  ce_survival <- function(u, cb) {
    i <- cbind(cb, seq_len(nd))
    mod$curve(L[i], U[i], pmax(u, K[i] * log10(1e-12)))
  }

  S <- matrix(NA_real_, n, nd); D <- matrix(0, n, nd)
  # Survival starts at the top of the curve at the first temperature where the
  # trace accrues injury, read as the constant shape reads its curve at zero
  # dose; a temperature at or below T_c, or a repeated time stamp, accrues none
  # and so does not set it.
  b1  <- bin[inj[1]]
  S[1, ] <- mod$curve(L[b1, ], U[b1, ], K[b1, ] * log10(1e-12))
  U1  <- U[b1, ]                             # repair scales with survival / U1 (at most 1)
  cap <- S[1, ]                              # lowest top met so far: deaths on exposure
  guess  <- rep(NA_real_, nd)
  u_cur <- rep(-Inf, nd); cb <- rep(b1, nd)   # cumulative exposure: state and its curve
  for (j in seq_len(n)[-1]) {
    a <- j - 1L; w <- dt[a]; b <- bin[a]
    if (w <= 0) {                            # a repeated time stamp: nothing changes
      S[j, ] <- S[a, ]; D[j, ] <- D[a, ]
      next
    }
    damaging <- !below[a]
    dD <- if (damaging) w / (10^MID[b, ] * unit_h) else 0   # relative dose this step
    if (damaging) {
      A[b, ] <- A[b, ] + dD
      if (!(b %in% rows)) rows <- c(rows, b)
      lo  <- pmin(lo, L[b, ]); hi <- pmin(hi, U[b, ])
      cap <- pmin(cap, U[b, ])
    }
    R <- 0                                   # repair this step, in relative doses
    if (!is.null(rep_rate)) {
      R <- rep_rate[a] * w * conv
      if (scale_with_surv) R <- R * S[a, ] / U1
      tot  <- colSums(A[rows, , drop = FALSE])
      frac <- ifelse(tot > 0, pmin(R / tot, 1), 0)
      A <- A * rep(1 - frac, each = nb)
      cleared <- frac >= 1                   # fully repaired: nothing outstanding
      A[, cleared] <- 0
      lo[cleared] <- Inf; hi[cleared] <- Inf
    }
    if (method == "individual") {
      root  <- rank_root(guess)
      guess <- root
      S[j, ] <- pmin(if (irreversible) S[a, ] else cap, root)
    } else {
      # A step that injures carries each draw's state to its curve (unless the
      # state is at or below that curve's floor) and adds the step's dose; at
      # or below T_c the state stays where it is. Repair lowers the dose on the
      # curve the state lies on. With irreversible mortality survival is the
      # lowest survival of the injury so far, so repaired injury kills no one
      # until it passes its previous peak.
      on_b <- rep(FALSE, nd)
      if (damaging) {
        mv <- which(cb != b)
        if (length(mv)) {
          m  <- ce_map(u_cur[mv], cb[mv], b, mv)
          go <- !m$stay
          u_cur[mv[go]] <- m$u[go]; cb[mv[go]] <- b
        }
        on_b <- cb == b
      }
      kc <- K[cbind(cb, seq_len(nd))]
      Dc <- pmax(0, 10^(u_cur / kc) + ifelse(on_b, dD, 0) - R)
      u_cur <- ifelse(Dc > 0, kc * log10(Dc), -Inf)
      S[j, ] <- pmin(if (irreversible) S[a, ] else cap, ce_survival(u_cur, cb))
    }
    D[j, ] <- if (is.null(p_target)) colSums(A[rows, , drop = FALSE]) else rank_injury(p_target)
  }
  list(survival = S, dose = D, max_iter = max_iter)
}

#' Predict heat injury and survival under a fluctuating temperature trace
#'
#' Propagates the model's full posterior through an Eulerian damage-
#' accumulation integration along the supplied temperature trace, returning
#' the posterior median and 95% credible band of:
#'
#' - **HI(t)** — cumulative heat injury, in percent of an LT_target_surv dose.
#'   When `HI(t) = 100`, the population has accumulated one full dose at the
#'   chosen survival threshold (by default midway between the asymptotes,
#'   `(low + up)/2`).
#' - **S(t)** — predicted survival fraction. With `shape = "constant"` it is
#'   read off the fitted dose-response curve (the log-logistic for
#'   [fit_4pl()], the log-normal or Weibull curve for [fit_lognormal()] /
#'   [fit_weibull()]) at the cumulative dose; with `shape = "varying"` it
#'   follows each individual's own heat injury (see the section "Varying
#'   shape").
#'
#' Optionally adds a temperature-dependent repair rate via
#' [repair_rate_schoolfield()].
#'
#' @param trace        Tibble with columns `time` (numeric time from start,
#'                     in `trace_unit`) and `temp` (°C), in time order.
#'                     Requires >= 2 rows.
#' @param workflow     Fitted `bayes_tls`.
#' @param target_surv  Threshold defining "1 dose". `"relative"` (default;
#'                     `(low + up)/2`), `"absolute"` (= 0.5), or a numeric in
#'                     `(0, 1)`. The default coincides with the classical LT50
#'                     when `low` is about 0 and `up` about 1; with sub-unit asymptotes
#'                     the relative threshold is the more biologically
#'                     meaningful anchor for dose accounting. With the
#'                     constant shape, a draw whose curve at the centring
#'                     temperature cannot reach an absolute target gives no
#'                     prediction (a warning reports how many). With
#'                     `shape = "varying"` and no repair, the target sets the
#'                     heat injury only; survival does not depend on it (with
#'                     repair it does, through the units of `r_ref`; see
#'                     `repair_pars`). In the varying shape (either
#'                     `varying_method`), heat injury for an absolute target is
#'                     the injury of the individual of that survival rank.
#'                     Where the target is at or above a trace temperature's
#'                     upper asymptote, that individual dies on exposure there
#'                     (above `T_c`), and heat injury is `Inf` (with a warning)
#'                     while injury accrued there is outstanding (repair that
#'                     clears all injury resets it); where it is at or below
#'                     the lower asymptote, that temperature adds no injury to
#'                     it.
#' @param T_c          Optional damage-accumulation threshold (°C). When
#'                     supplied, the damage rate is forced to zero at
#'                     `temp <= T_c` (matches the heat-injury integral,
#'                     Equation 8 of the manuscript).
#'                     Default `NULL` lets the rate fall naturally with T.
#' @param trace_unit   Time unit of the trace's `time` column: one of
#'                     `"hours"` (default), `"minutes"`, `"seconds"`, `"days"`.
#'                     Reconciled internally with the model's fitted
#'                     `duration_unit`, so the result is correct for any
#'                     combination of model and trace time units.
#' @param ndraws       Posterior draws to use; `NULL` (default) uses the full
#'                     posterior, keeping the integral consistent with the model.
#'                     Pass an integer to subsample -- useful for speed, since the
#'                     full posterior over a long temperature trace can be slow.
#' @param repair       Logical. If `TRUE`, add Sharpe-Schoolfield repair.
#'                     Default `FALSE`.
#' @param repair_pars  Required when `repair = TRUE`. Named list with elements
#'                     `TA, TAL, TAH, TL, TH, TREF, r_ref` passed straight to
#'                     [repair_rate_schoolfield()]. `r_ref` is in doses per
#'                     hour, the damage-rate units, where a dose is the
#'                     threshold time of `target_surv`; a given `r_ref`
#'                     therefore repairs more or less depending on the target.
#'                     (With `shape = "varying"`, injury is kept in doses of
#'                     the relative threshold, `r_ref` is converted with the
#'                     two thresholds' ratio at the centring temperature, as in
#'                     the constant shape; repair removes the same fraction of
#'                     every individual's injury, or with
#'                     `varying_method = "cumulative"` lowers the shared
#'                     injury by the repair amount.)
#' @param repair_scales_with_survival Logical. If `TRUE` (default), repair
#'                     rate at each step is scaled by `survival / up` so dead
#'                     organisms do not contribute to repair. (With
#'                     `shape = "varying"`, `up` is that of the first
#'                     temperature that accrues injury, where survival starts,
#'                     so the scale is at most 1.)
#' @param irreversible_mortality Logical. If `TRUE` (default), survival can
#'                     only decrease over time — i.e. once the population's
#'                     predicted survival reaches a value, it cannot rebound
#'                     even if cumulative dose subsequently decreases.
#' @param save_draws   Logical. If `TRUE`, return the full per-draw
#'                     trajectories. Default `FALSE`.
#' @param seed         Optional integer seeding the posterior-draw subsample for
#'                     reproducibility. `NULL` (default) leaves the RNG untouched.
#' @param by           Optional moderator column(s) for per-group injury. `NULL`
#'                     (default) uses the fit's moderators; a single-condition fit
#'                     returns the ungrouped result. A grouped fit runs the dose
#'                     integral through each group's 4PL and `summary` gains the
#'                     moderator column(s).
#' @param shape        How the curve's shape enters the prediction.
#'                     `"constant"` (default) holds `low`, `up` and `k` at the
#'                     centring temperature `T_bar` and lets only the midpoint
#'                     shift with temperature, along a straight line through
#'                     its value at `T_bar`: survival is read off that one
#'                     curve at the accumulated dose, the classical
#'                     constant-shape heat-injury assumption (via
#'                     [hi_pars()]). `"varying"` (experimental) uses the
#'                     temperature-local `low(T)`, `up(T)`, `k(T)` and
#'                     `mid(T)` (via [extract_4pl_pars()]); `varying_method`
#'                     sets how survival follows them. See the section
#'                     "Varying shape".
#' @param varying_method For `shape = "varying"` only: `"individual"`
#'                     (default; each individual accumulates injury against
#'                     its own threshold time) or `"cumulative"` (cumulative
#'                     exposure: survival is the state). See the section
#'                     "Varying shape".
#' @return A list with elements:
#'   - `summary`: tibble with `time`, `temp`, and posterior median + 95%
#'     CrI for `hi`, `survival`, and `mortality` at each time step (plus the
#'     moderator column(s) for a grouped fit).
#'   - `draws`: optional per-draw trajectories (when `save_draws = TRUE`); a
#'     draw that gives no prediction is kept as rows of `NA`.
#'   - `meta`: list of inputs used. `meta$ndraws` is the number of posterior
#'     draws that produced a prediction; a warning reports any that did not,
#'     and they are left out of every summary column. For a grouped fit it is
#'     the smallest such number in any group, and `meta$ndraws_by_group`
#'     gives each group's count.
#' @section Varying shape (experimental):
#' The fitted curves describe survival at constant temperatures. When their
#' shape (`low`, `up`, `k`) changes with temperature, survival under a
#' fluctuating trace depends on why it changes, which constant-temperature data
#' cannot reveal. `varying_method` chooses between two coherent methods. Both
#' reproduce the fitted curve at a constant temperature; both use the constant
#' shape's rule that, with irreversible mortality, injury that is repaired and
#' re-accrued kills no one until it passes its previous peak; and both equal the
#' constant shape exactly when `low`, `up` and `k` do not change with
#' temperature and the midpoint is linear in it (with an absolute target and a
#' very shallow curve, `k` below about 2, the starting values differ slightly:
#' the two shapes read their curves at a tiny dose measured in different
#' units). Otherwise they differ, and comparing them shows how much the unknown
#' mechanism matters.
#'
#' - `"individual"` (default): each individual is followed separately. The
#'   individual of survival rank `p` is the one that dies, at a constant
#'   temperature, when survival there falls to `p`, at that temperature's
#'   threshold time for `p`. It accumulates heat injury additively (the time at
#'   each temperature divided by its own threshold time there) and dies when
#'   its injury reaches 1. This is exact if individuals have their own TDT lines
#'   that fan out and keep their ranks. Without repair, once the trace has begun
#'   to accrue injury, survival depends only on the time spent at each
#'   temperature, not on their order or on where the trace starts; more time at
#'   any temperature never raises it (a hotter trace can still end higher where
#'   the fitted curves cross); and the heat injury of an absolute target reaches
#'   100% exactly when survival falls to that target. Repair removes the same
#'   fraction of every individual's injury.
#' - `"cumulative"`: cumulative exposure, with survival as the state. Over each
#'   time step, survival continues along the current temperature's curve from
#'   the point where that curve has the current survival. This is exact if
#'   individuals share one injury process whose rate depends on the injury
#'   already accrued, and differ only in their thresholds. It is the usual rule
#'   in reliability for a shape that changes with stress (Nelson 1980) and the
#'   one used by Rezende et al. (2020) and Bullard et al. (2026); given the same
#'   curves it reproduces Bullard et al.'s code, except that their code reads
#'   each step's curve at the step's end rather than its start (the difference
#'   shrinks as the steps do). Temperature order matters. Because survival is
#'   the state, deaths on exposure to a lower `up` (below) also advance the
#'   survivors along later curves. A step injures only while survival is above
#'   the floor of its curve; otherwise the state stays on the curve of the last
#'   step that injured it (as ranks below `low` accrue no injury in the
#'   per-individual method), as it does at or below `T_c`. With repair, a
#'   shared injury falls by the repair amount on the dose scale of the curve the
#'   state lies on, so the repair during a cool spell depends on the curve of
#'   the last step that injured. The heat injury of an absolute target is the
#'   per-individual method's and does not mark when survival under this method
#'   reaches the target.
#'
#' Each time step is held at the temperature at its start, as the dose is.
#' Survival starts at `up` of the first temperature at which the trace accrues
#' injury (above `T_c`); a trace that accrues none keeps the constant shape's
#' value, `up` at `T_bar`. Both
#' methods read a lower `up` at some temperature as individuals that die as
#' soon as they are exposed to it (a convention without published precedent),
#' so survival drops to that `up` once the trace spends time there above
#' `T_c`; those deaths are permanent even with
#' `irreversible_mortality = FALSE`, as the share `1 - up` is in the constant
#' shape. If a temperature-dependent `up` reflects losses unrelated to heat, fit
#' `up` as a constant or set `T_c`; `up(T)` is extrapolated at trace
#' temperatures outside the assay range. Heat injury is the same under both
#' methods unless repair scales with survival, when each method's own survival
#' scales it. For a continuous response (the Beta family) the curve is a mean
#' response, so the reading in terms of individuals does not apply. With
#' `irreversible_mortality = FALSE`, survival rebounds as injury is repaired,
#' up to the lowest `up` met so far. The cost of the varying shape grows with
#' the number of distinct trace temperatures, so rounding the trace (for
#' example to 0.1 °C) speeds it up.
#' @references
#' Bullard GW, Buckley LB, Kingsolver JG (2026). Incorporating variation in
#' death times improves predictions of ectotherm responses to stressful
#' temperatures. PLOS Biology 24(5): e3003623. \doi{10.1371/journal.pbio.3003623}
#'
#' Nelson W (1980). Accelerated life testing - step-stress models and data
#' analyses. IEEE Transactions on Reliability R-29(2): 103-108.
#' \doi{10.1109/TR.1980.5220742}
#'
#' Rezende EL, Bozinovic F, Szilagyi A, Santos M (2020). Predicting temperature
#' mortality and selection in natural Drosophila populations. Science
#' 369(6508): 1242-1245. \doi{10.1126/science.aba9287}
#' @examples
#' \dontrun{
#' scens <- make_temperature_scenarios()
#' hi    <- predict_heat_injury(scens$single_spike, wf)
#' hi$summary
#' }
#' @export
predict_heat_injury <- function(trace, workflow,
                                target_surv = "relative",
                                T_c         = NULL,
                                trace_unit  = "hours",
                                ndraws      = NULL,
                                repair      = FALSE,
                                repair_pars = NULL,
                                repair_scales_with_survival = TRUE,
                                irreversible_mortality      = TRUE,
                                save_draws                  = FALSE,
                                seed                        = NULL,
                                by                          = NULL,
                                shape = c("constant", "varying"),
                                varying_method = c("individual", "cumulative")) {
  shape <- match.arg(shape)
  varying_method <- match.arg(varying_method)
  if (!has_fit(workflow))
    stop("workflow$fit is NULL. Fit the model first.", call. = FALSE)
  tdt_check_ndraws(ndraws)
  if (nrow(trace) < 2L)
    stop("trace must have at least 2 rows.", call. = FALSE)
  if (!all(c("time", "temp") %in% names(trace)))
    stop("trace must have `time` and `temp` columns.", call. = FALSE)
  if (anyNA(trace$time) || anyNA(trace$temp))
    stop("trace has NA in `time` or `temp`. Interpolate or drop missing rows ",
         "first -- an NA temperature would otherwise be silently counted as zero ",
         "heat injury, under-counting the dose.", call. = FALSE)
  if (repair && is.null(repair_pars))
    stop("Supply repair_pars when repair = TRUE.", call. = FALSE)

  # Normalise the threshold up front so we can pass either the literal
  # probability or the sentinel string "relative" down to the helpers.
  ts <- resolve_target_surv(target_surv)
  ts_arg <- if (ts$mode == "relative") "relative" else ts$prob

  trace     <- trace[order(trace$time), , drop = FALSE]
  n         <- nrow(trace)
  temp_mean <- workflow$meta$temp_mean
  # The fitted dose-response model, resolved once: the integrators below call
  # its threshold and dose-to-survival maps for every draw and time step.
  mod       <- tdt_model_spec(tdt_model(workflow))

  # Reconcile BOTH time units to a common base (hours), so the integral is
  # correct no matter how the model was fit OR how the trace is labelled:
  #   - `tau` is returned in the model's `duration_unit`;
  #   - the trace step is in `trace_unit` (the `time` column).
  # We convert the damage rate to "doses per hour" and the step `dt` to hours.
  # Without this, e.g. a minutes-fitted model driven by an hours trace would
  # under-count the accumulated dose 60-fold (and the reverse over-counts).
  to_hours <- function(u) switch(as.character(u),
                   seconds = 1 / 3600, minutes = 1 / 60, hours = 1, days = 24,
                   stop(sprintf(
                     "predict_heat_injury(): unsupported time unit '%s'; expected one of seconds/minutes/hours/days.",
                     u), call. = FALSE))
  unit_h <- to_hours(workflow$meta$duration_unit %||% "hours")               # model -> hours
  # Per-interval widths (hours). Using diff() rather than reusing a single first
  # step is what makes irregular / gappy traces integrate correctly; dt_vec[j-1]
  # is the width of the interval ending at point j.
  dt_vec <- diff(trace$time) * to_hours(trace_unit)                          # trace -> hours

  by       <- tdt_resolve_by(workflow, by)
  pars_all <- hi_pars(workflow, by = by)
  # ndraws = NULL (default) uses the FULL posterior -- keeps the integral
  # consistent with the model and the other derivation functions; pass an
  # integer to subsample for speed on long traces.
  if (is.null(ndraws)) ndraws <- nrow(pars_all)

  # Per-draw Euler dose-accumulation integral for one set of 4PL parameter draws.
  integrate_pars <- function(pars) {
    pred_list <- vector("list", nrow(pars))
    for (i in seq_len(nrow(pars))) {
      tau <- time_to_surv_threshold_4pl(
        temp = trace$temp, survival_target = ts_arg,
        low = pars$low[i], up = pars$up[i], k = pars$k[i],
        mid_int = pars$mid_int[i], mid_temp = pars$mid_temp[i],
        temp_mean = temp_mean, model = mod
      )
      dmg <- 1 / (tau * unit_h)               # doses per hour (unit-reconciled)
      dmg[!is.finite(dmg)] <- 0
      if (!is.null(T_c)) dmg[trace$temp <= T_c] <- 0

      rep_rate <- if (repair) {
        repair_rate_schoolfield(
          temp_celsius = trace$temp,
          TA = repair_pars$TA, TAL = repair_pars$TAL,
          TAH = repair_pars$TAH, TL = repair_pars$TL,
          TH = repair_pars$TH, TREF = repair_pars$TREF,
          r_ref = repair_pars$r_ref
        )
      } else {
        rep(0, n)
      }

      # Forward-Euler dose ODE: d(dose)/dt = damage(T) - repair(T) * scale. The
      # dose at the first time point is zero (no exposure has elapsed); each later
      # point adds the rate at the START of its interval times that interval's
      # width. (The previous loop credited a full step at t = 0 and reused the
      # first dt for every interval -- over-counting by one step and corrupting
      # irregular traces.) Kept consistent with planted_dose_from_trace().
      dose <- numeric(n); survival <- numeric(n)
      dose[1]     <- 0
      survival[1] <- survival_from_dose(0, low = pars$low[i], up = pars$up[i],
                                        k = pars$k[i], target_surv = ts_arg,
                                        model = mod)
      for (j in seq_len(n)[-1]) {
        w     <- dt_vec[j - 1]
        rep_j <- rep_rate[j - 1] * w
        if (repair_scales_with_survival) rep_j <- rep_j * survival[j - 1] / pars$up[i]
        new_dose <- max(0, dose[j - 1] + dmg[j - 1] * w - rep_j)
        surv_raw <- survival_from_dose(
          new_dose, low = pars$low[i], up = pars$up[i], k = pars$k[i],
          target_surv = ts_arg, model = mod
        )
        survival[j] <- if (irreversible_mortality) min(survival[j - 1], surv_raw) else surv_raw
        dose[j]     <- new_dose
      }
      pred_list[[i]] <- data.frame(
        .draw = pars$.draw[i], time = trace$time, temp = trace$temp,
        dose = dose, hi = dose * 100, survival = survival, mortality = 1 - survival
      )
    }
    dplyr::bind_rows(pred_list)
  }

  # ----- shape = "varying": temperature-LOCAL low/up/k/mid -----
  # The local curves are evaluated once per group at the trace's unique temps;
  # `curve$.draw` indexes the full posterior, so a sampled draw is matched by its
  # .draw id. `tcol` maps each trace point to its temperature column.
  temps_unique <- sort(unique(trace$temp))
  tcol         <- match(trace$temp, temps_unique)
  curves_by_grp <- if (shape == "varying")
    hi_local_curves(workflow, by = by, temps_unique = temps_unique) else NULL
  relative_target <- is.character(ts_arg) && length(ts_arg) == 1L &&
    ts_arg == "relative"
  below_tc <- if (is.null(T_c)) rep(FALSE, n) else trace$temp <= T_c
  # Problems met while integrating the varying shape, reported once at the end:
  # draws whose absolute target is already passed at some trace step, the
  # temperatures where that happens, draws whose local curve is invalid at a
  # trace temperature, and draws whose repair rate cannot be converted.
  vary_issues <- list(passed = 0L, passed_temps = numeric(0), invalid = 0L,
                      no_conv = 0L)

  # Varying-shape integrator: individual thresholds (see hi_rank_integrate()).
  # All draws of a group are integrated together. A draw with no prediction
  # keeps its rows, as NA, so its group still appears in the summary and the
  # draw is counted as dropped.
  integrate_pars_varying <- function(pars, curve) {
    nd   <- nrow(pars)
    rmap <- match(pars$.draw, curve$.draw)    # row of each draw in the local curves
    keep <- which(!is.na(rmap))
    vary_issues$invalid <<- vary_issues$invalid + (nd - length(keep))
    rep_rate <- if (repair) {
      repair_rate_schoolfield(
        temp_celsius = trace$temp,
        TA = repair_pars$TA, TAL = repair_pars$TAL,
        TAH = repair_pars$TAH, TL = repair_pars$TL,
        TH = repair_pars$TH, TREF = repair_pars$TREF,
        r_ref = repair_pars$r_ref)
    } else NULL
    # `r_ref` is in doses of the target per hour; injury is kept in relative
    # doses, so convert with the target's offset at the centring temperature,
    # as the constant shape implies.
    conv <- rep(1, nd)
    if (!relative_target) {
      conv <- 10^tdt_abs_offset(mod, pars$low, pars$up, pars$k, ts_arg)
      if (repair) {
        no_conv <- keep[!is.finite(conv[keep])]
        vary_issues$no_conv <<- vary_issues$no_conv + length(no_conv)
        keep <- setdiff(keep, no_conv)
      }
    }
    surv <- matrix(NA_real_, n, nd); dose <- matrix(NA_real_, n, nd)
    if (length(keep)) {
      r <- rmap[keep]
      L <- t(curve$low[r, , drop = FALSE]); U <- t(curve$up[r, , drop = FALSE])
      K <- t(curve$k[r, , drop = FALSE]);   M <- t(curve$mid[r, , drop = FALSE])
      res <- hi_rank_integrate(
        L, U, K, M, bin = tcol, dt = dt_vec, below = below_tc, unit_h = unit_h,
        mod = mod, p_target = if (relative_target) NULL else ts_arg,
        rep_rate = rep_rate, conv = conv[keep],
        scale_with_surv = repair_scales_with_survival,
        irreversible = irreversible_mortality, method = varying_method,
        s_none = mod$curve(pars$low[keep], pars$up[keep], pars$k[keep] * log10(1e-12)))
      surv[, keep] <- res$survival; dose[, keep] <- res$dose
      if (!relative_target) {
        # Steps that take time, above T_c, where the target is at or above up.
        used <- which(dt_vec > 0 & !below_tc[-n])
        if (length(used)) {
          passed <- ts_arg >= U[tcol[used], , drop = FALSE]
          if (any(passed)) {
            vary_issues$passed <<- vary_issues$passed + sum(colSums(passed) > 0)
            vary_issues$passed_temps <<- unique(c(vary_issues$passed_temps,
                                                  trace$temp[used][rowSums(passed) > 0]))
          }
        }
      }
    }
    data.frame(.draw = rep(pars$.draw, each = n), time = rep(trace$time, nd),
               temp = rep(trace$temp, nd), dose = as.vector(dose),
               hi = 100 * as.vector(dose), survival = as.vector(surv),
               mortality = 1 - as.vector(surv))
  }

  # Dispatch: pick the per-draw integrator for the requested shape. The grouped
  # path selects the matching group's local curves by label.
  run_integral <- function(pars, grp_label) {
    if (shape == "constant") integrate_pars(pars)
    else integrate_pars_varying(pars, curves_by_grp[[grp_label]])
  }
  # Draws that produced a complete survival trajectory. A draw gives no
  # prediction when an absolute target lies outside its curve at the centring
  # temperature (constant shape), or, for the varying shape, when its local
  # curve is invalid at a trace temperature or its repair rate cannot be
  # converted. `meta$ndraws` reports these draws, not the number requested.
  count_used <- function(d) {
    if (is.null(d) || !nrow(d)) return(0L)
    as.integer(sum(tapply(is.finite(d$survival), d$.draw, all)))
  }
  # A draw with no complete survival trajectory is blanked in every column,
  # so it is left out of the heat-injury summary as well as the survival one.
  blank_dropped <- function(d) {
    if (is.null(d) || !nrow(d)) return(d)
    ok  <- tapply(is.finite(d$survival), d$.draw, all)
    bad <- as.character(d$.draw) %in% names(ok)[!ok]
    d[bad, c("dose", "hi", "survival", "mortality")] <- NA_real_
    d
  }

  if (!is.null(seed)) set.seed(seed)   # reproducible posterior-draw subsample
  if (is.null(by)) {
    n_req <- min(ndraws, nrow(pars_all))
    draws <- blank_dropped(run_integral(dplyr::slice_sample(pars_all, n = n_req), "all"))
    n_used     <- count_used(draws)
    n_by_group <- NULL
    n_req_all  <- n_req
    n_used_all <- n_used
    empty_groups <- character(0)
  } else {
    groups <- unique(pars_all[, by, drop = FALSE])
    per_group_n    <- integer(nrow(groups))
    per_group_used <- integer(nrow(groups))
    grp_labels     <- character(nrow(groups))
    draws <- dplyr::bind_rows(lapply(seq_len(nrow(groups)), function(gi) {
      gp     <- dplyr::inner_join(pars_all, groups[gi, , drop = FALSE], by = by)
      n_gi   <- min(ndraws, nrow(gp))
      per_group_n[gi] <<- n_gi
      # Group label, built EXACTLY as tls_build_grid() builds its `.grp`
      # (do.call(paste, c(<by columns>, sep = " / "))) so the varying-mode local
      # curves -- keyed by that same `.grp` -- are matched for this group.
      grp_label <- do.call(paste, c(groups[gi, by, drop = FALSE], sep = " / "))
      grp_labels[gi] <<- grp_label
      res <- blank_dropped(run_integral(dplyr::slice_sample(gp, n = n_gi), grp_label))
      per_group_used[gi] <<- count_used(res)
      if (is.null(res) || !nrow(res)) return(NULL)
      cbind(groups[gi, , drop = FALSE], res, row.names = NULL)
    }))
    # The smallest number of draws behind any group's summary; each group's
    # count is in `n_by_group`.
    n_by_group <- stats::setNames(per_group_used, grp_labels)
    n_used     <- min(per_group_used)
    n_req_all  <- sum(per_group_n)
    n_used_all <- sum(per_group_used)
    empty_groups <- grp_labels[per_group_used == 0L]
  }

  if (n_used_all < n_req_all) {
    why <- if (shape == "constant") {
      "the survival target lies outside those draws' fitted curve at the centring temperature"
    } else {
      paste(c(if (vary_issues$invalid > 0)
                sprintf("%d had an invalid local curve at a trace temperature",
                        vary_issues$invalid),
              if (vary_issues$no_conv > 0)
                sprintf(paste("%d had the survival target outside their curve at the",
                              "centring temperature, so the repair rate (in doses of",
                              "the target) could not be converted"), vary_issues$no_conv)),
            collapse = "; ")
    }
    warning(sprintf(paste("predict_heat_injury(): %d of %d posterior draws gave no",
                          "prediction and are left out of the summary (%s).%s",
                          "`meta$ndraws` counts the draws used."),
                    n_req_all - n_used_all, n_req_all, why,
                    if (length(empty_groups))
                      sprintf(" No draw gave a prediction for %s, so its summary is NA.",
                              paste(empty_groups, collapse = ", "))
                    else ""), call. = FALSE)
  }
  if (vary_issues$passed > 0) {
    tt <- sort(vary_issues$passed_temps)
    tt_txt <- if (length(tt) <= 4L) paste(format(tt), collapse = ", ")
              else sprintf("%s to %s", format(min(tt)), format(max(tt)))
    warning(sprintf(paste("predict_heat_injury(): the survival target %s is at or",
                          "above the fitted curve's upper asymptote at %s \u00b0C in",
                          "%d of %d draws, so it is already passed when exposure",
                          "starts there. Heat injury relative to this target is",
                          "infinite (Inf) while injury accrued there is outstanding;",
                          "the survival prediction is not affected. Use",
                          "target_surv = \"relative\" or a lower target for finite",
                          "heat injury."),
                    format(ts$prob), tt_txt, vary_issues$passed, n_req_all),
            call. = FALSE)
  }

  q_lower <- function(x) stats::quantile(x, 0.025, na.rm = TRUE)
  q_upper <- function(x) stats::quantile(x, 0.975, na.rm = TRUE)
  summary <- draws |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(by, "time", "temp")))) |>
    dplyr::summarise(
      hi_median   = stats::median(hi,        na.rm = TRUE),
      hi_lower    = q_lower(hi),
      hi_upper    = q_upper(hi),
      surv_median = stats::median(survival,  na.rm = TRUE),
      surv_lower  = q_lower(survival),
      surv_upper  = q_upper(survival),
      mort_median = stats::median(mortality, na.rm = TRUE),
      mort_lower  = q_lower(mortality),
      mort_upper  = q_upper(mortality),
      .groups     = "drop"
    )

  out <- list(
    summary = summary,
    meta    = list(target_surv = ts$label, target_mode = ts$mode,
                   target_prob = ts$prob, T_c = T_c, repair = repair,
                   repair_pars = repair_pars, ndraws = n_used,
                   ndraws_by_group = n_by_group, by = by, shape = shape,
                   varying_method = if (shape == "varying") varying_method)
  )
  if (save_draws) out$draws <- draws
  out
}
