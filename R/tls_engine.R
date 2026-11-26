# Shared TDT extraction engine. Every post-fit quantity (z, CTmax, T_crit, LT
# curves, survival, heat injury) is derived from brms::posterior_linpred(nlpar=)
# evaluated at an explicit grid — NEVER from raw b_*_Intercept coefficient names.
# posterior_linpred(nlpar = "mid") returns the linear predictor for a midpoint
# fit and the nlf-derived mid for a direct fit, so this engine is identically
# parameterisation- and coding-agnostic. tls() and the single-condition readers
# (derive_z, derive_temperature_for_duration, tdt_parameter_table,
# extract_4pl_pars) all build on these helpers; the pure-math layer
# (tls_local_z / tls_invert_logLT) operates on already-evaluated matrices so the
# z/CTmax math stays fast-testable without a Stan fit.

# Build a moderator x temperature grid for posterior_linpred. Crosses `by`-levels
# with `temp_grid`, fills any model-data columns the grid lacks (so brms
# validates), and tags groups in `.grp` (`"all"` when `by = NULL`). `temp` is the
# (centred) temperature column name in the model data.
tls_build_grid <- function(mdata, by = NULL, temp = "temp_c", temp_grid = NULL,
                           newdata = NULL, censored = FALSE) {
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
  # A censored fit's flag only changes how a reading enters the likelihood; the
  # grid reads the curve, so mark its rows uncensored rather than copying row 1.
  # Only for a censored fit (`censored = TRUE`), and never a moderator column:
  # an uncensored fit may carry its own column of that name.
  if (isTRUE(censored) && "censored" %in% names(newdata) && !"censored" %in% by)
    newdata$censored <- "none"
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

# Validate a user-facing `ndraws`: NULL (all draws) or one whole number >= 1.
# brms used to check this when `ndraws` was passed straight to
# posterior_linpred(); the readers now subsample through tls_draw_ids().
tdt_check_ndraws <- function(ndraws) {
  if (is.null(ndraws)) return(invisible(NULL))
  if (!is.numeric(ndraws) || length(ndraws) != 1L || !is.finite(ndraws) ||
      ndraws < 1 || ndraws != round(ndraws))
    stop("`ndraws` must be NULL (all draws) or a single whole number >= 1.",
         call. = FALSE)
  invisible(ndraws)
}

# Drop non-finite temperatures from a user temp_grid, with a warning; error if
# none remain. NULL passes through (the caller's default grid). An Inf left in
# a grid would otherwise make every crossing search fail silently.
tdt_finite_grid <- function(temp_grid) {
  if (is.null(temp_grid)) return(NULL)
  fin <- is.finite(temp_grid)
  if (!any(fin)) stop("`temp_grid` has no finite temperatures.", call. = FALSE)
  if (!all(fin))
    warning("Dropping ", sum(!fin), " non-finite value(s) from `temp_grid`.",
            call. = FALSE)
  temp_grid[fin]
}

# Subsample draw indices reproducibly, matching the retired tdt_extract_pars()
# selection (`d[sort(sample.int(n, ndraws)), ]`) so the single-condition readers
# stay bit-identical to the coefficient-parsing implementation. Returns NULL
# (use all draws, in order) when `ndraws` is NULL or >= the available draws. The
# caller is responsible for set.seed() beforehand.
tls_draw_ids <- function(fit, ndraws) {
  tdt_check_ndraws(ndraws)
  if (is.null(ndraws)) return(NULL)
  n <- brms::ndraws(fit)
  if (ndraws >= n) return(NULL)
  sort(sample.int(n, ndraws))
}

# Asymptote windows of the disjoint-bounds reparameterisation for a workflow:
# the ones it was fitted with (meta$bounds), else rebuilt from its response
# bounds (meta$lower / meta$upper, default 0 / 1), as tls() does. A workflow
# from fit_4pl() always carries meta$bounds; a hand-built bayes_tls may not.
tdt_workflow_bounds <- function(workflow) {
  workflow$meta$bounds %||%
    compute_4pl_bounds(workflow$meta$lower %||% 0, workflow$meta$upper %||% 1)
}

# Evaluate the four 4PL sub-parameters at `newdata`, returning per-draw matrices
# (`[ndraws x nrow(newdata)]`) on the natural scale plus the threshold-specific
# logLT (relative: mid; absolute p: mid + the model's offset, see
# tdt_abs_offset()). This is the
# ONLY engine layer that calls brms — the gate boundary for fast vs gated tests.
# `draw_ids` selects specific posterior draws (see tls_draw_ids), shared by all
# four posterior_linpred() calls. There is deliberately no `ndraws` argument:
# brms draws a fresh random subsample on every call, which would pair the four
# sub-parameters from different posterior draws. `warn = FALSE` silences the
# out-of-asymptote warning (for tls_extend_crossings(), whose extrapolated
# cells report through the caller's own warning).
tls_eval_subpars <- function(fit, newdata, bounds,
                             nlpars = c("lowraw", "upraw", "logk", "mid"),
                             re_formula = NA, draw_ids = NULL,
                             mode = "absolute", p = 0.5, warn = TRUE) {
  need <- c("low_min", "low_w", "up_min", "up_w")
  if (!is.list(bounds) ||
      !all(vapply(need, function(nm) is.numeric(bounds[[nm]]) &&
                    length(bounds[[nm]]) == 1L, logical(1))))
    stop("Internal error: `bounds` must be a compute_4pl_bounds() list ",
         "(low_min, low_w, up_min, up_w). For a workflow use ",
         "tdt_workflow_bounds(workflow).", call. = FALSE)
  lp <- function(np) brms::posterior_linpred(fit, newdata = newdata, nlpar = np,
                                             re_formula = re_formula,
                                             draw_ids = draw_ids)
  low <- bounds$low_min + stats::plogis(lp(nlpars[1])) * bounds$low_w
  up  <- bounds$up_min  + stats::plogis(lp(nlpars[2])) * bounds$up_w
  k   <- exp(lp(nlpars[3]))
  mid <- lp(nlpars[4])
  if (identical(mode, "absolute")) {
    # Absolute LT at survival p: the midpoint plus the fitted model's offset
    # (tdt_abs_offset(); for the log-logistic log((up - p)/(p - low))/k, as
    # before). The model is read from the fit itself, so a bare brmsfit and
    # every caller -- tls(), derive_z(), derive_temperature_for_duration() --
    # get the right curve. p is attainable only where the fitted curve spans
    # it, low < p < up; elsewhere the LT at p does not exist for that
    # draw x temperature cell, so the cell is NA (no bare "NaNs produced") and
    # we warn ONCE with the reason.
    offset <- tdt_abs_offset(tdt_model(fit), low, up, k, p)
    bad    <- is.na(offset)
    logLT  <- mid + offset
    if (any(bad) && isTRUE(warn))
      warning(tls_abs_outside_message(bad, low, up, p), call. = FALSE)
  } else {
    logLT <- mid
  }
  list(low = low, up = up, k = k, mid = mid, logLT = logLT)
}

# Warning text for an absolute threshold p lying outside the fitted asymptotes
# (pure; raised by tls_eval_subpars()). `bad`, `low`, `up` are the
# [draws x evaluation points] matrices; the points can repeat a temperature
# (tls() evaluates each one three times), so the message counts draws and
# gives the share of evaluations, with the side p falls on.
tls_abs_outside_message <- function(bad, low, up, p) {
  bad <- as.matrix(bad)
  n_bad <- sum(bad)
  n_up  <- sum(bad & !is.na(up)  & up  <= p)
  n_low <- sum(bad & !is.na(low) & low >= p)
  n_na  <- n_bad - n_up - n_low
  sides <- sprintf(paste0("the upper asymptote is at or below %.3g (survival ",
                          "starts below the target) in %d of them and the ",
                          "lower asymptote at or above %.3g (survival never ",
                          "falls to it) in %d"), p, n_up, p, n_low)
  if (n_na > 0L) sides <- sprintf("%s; %d have undefined asymptotes", sides, n_na)
  sprintf(paste0(
    "target_surv = %.3g (absolute) is outside the fitted curve's asymptotes for ",
    "%d of %d draws at one or more of the temperatures evaluated (%d of %d ",
    "evaluations, %s%%): %s. The LT at %.3g survival does not exist there. Use a ",
    "threshold within the asymptote range, or target_surv = \"relative\" (the ",
    "per-draw midpoint, always defined)."),
    p, sum(rowSums(bad) > 0L), nrow(bad), n_bad, length(bad),
    formatC(100 * n_bad / length(bad), digits = 2, format = "fg"), sides, p)
}

# --- pure math on already-evaluated matrices (no brms; fast-testable) ---------

# Per draw, TRUE when the evaluated logLT is a straight line in temperature:
# every point lies on the line through the lowest- and highest-temperature
# points, to within `tol` of the line's rise across them (rounding error is
# ~1e-12 of it). The reference-exposure crossing then has an exact closed
# form at any temperature, beyond the evaluated grid too. This holds for every
# relative threshold of a 4PL whose midpoint is linear in temperature
# (midpoint and direct fits alike) and for an absolute threshold when the
# asymptotes and steepness do not change with temperature. NA for a draw with
# a non-finite value. `M` is `[ndraws x length(temps)]`; `temps` need not be
# sorted or unique.
tls_linear_draws <- function(M, temps, tol = 1e-7) {
  lo <- which.min(temps); hi <- which.max(temps)
  span <- temps[hi] - temps[lo]
  if (!is.finite(span) || span <= 0) return(rep(NA, nrow(M)))
  slope <- (M[, hi] - M[, lo]) / span
  pred  <- M[, lo] + outer(slope, temps - temps[lo])
  dev   <- apply(abs(M - pred), 1L, max)
  lin   <- dev <= tol * pmax(abs(slope) * span, 1e-6)
  lin[!is.finite(dev) | !is.finite(slope)] <- NA
  lin
}

# The temperatures a bent LT curve is searched over beyond the requested grid:
# the sorted grid `tg` plus points beyond each end, out to `width` (at least
# the grid's own width and 1 degree) at a quarter of the grid's median spacing
# (`refine`; at most `max_per_side` points a side, which then sets the
# spacing). Inside the grid the points -- and so any crossing found there --
# are unchanged.
tls_extend_grid <- function(tg, width = NULL, refine = 4L, max_per_side = 80L) {
  tg <- sort(unique(tg))
  w  <- max(width %||% 0, diff(range(tg)), 1)
  s  <- (if (length(tg) > 1L) stats::median(diff(tg)) else w) / refine
  n  <- min(max_per_side, max(1L, ceiling(w / s - 1e-9)))
  s  <- w / n
  c(min(tg) - s * rev(seq_len(n)), tg, max(tg) + s * seq_len(n))
}

# The grid a bent LT curve is inverted on: the sorted grid `tg` with points
# inserted in every interval so that consecutive points are at most `step`
# degrees apart (at most `max_refine` sub-intervals per interval). Linear
# interpolation between grid points then puts a crossing within a small
# fraction of the step that would otherwise bound its error (up to ~0.1-2 C
# per draw on a coarse grid next to an asymptote). The original points are
# kept exactly; a grid already that fine is returned unchanged.
tls_refine_grid <- function(tg, step = 0.05, max_refine = 20L) {
  tg <- sort(unique(tg))
  if (length(tg) < 2L) return(tg)
  out <- tg[1]
  for (i in seq_len(length(tg) - 1L)) {
    d <- tg[i + 1L] - tg[i]
    n <- min(max_refine, max(1L, ceiling(d / step - 1e-9)))
    out <- c(out, if (n > 1L) tg[i] + d * seq_len(n - 1L) / n, tg[i + 1L])
  }
  out
}

# Sign-safe logLT for inverting an absolute-threshold LT curve. Where the
# target p lies outside a draw's fitted asymptotes the LT at p does not exist:
# survival is below p at every duration when up <= p (an LT of zero) and above
# p at every duration when low >= p (an infinite LT). Those cells get logLT =
# -Inf / +Inf, so the sign of logLT - target is right everywhere; NA stays NA
# (undefined sub-parameters). `sp` is a tls_eval_subpars() list; a relative
# threshold (p = NA) returns sp$logLT unchanged.
tls_clamp_logLT <- function(sp, p) {
  L <- sp$logLT
  if (is.na(p)) return(L)
  L[is.na(L) & !is.na(sp$up)  & sp$up  <= p] <- -Inf
  L[is.na(L) & !is.na(sp$low) & sp$low >= p] <-  Inf
  L
}

# Survival at log10 duration `target`, minus p. It is finite wherever the
# sub-parameters are and has the sign of logLT - target (survival at t is above
# p exactly when the LT at p is longer than t), so it places a crossing next to
# a cell where logLT is infinite (tls_clamp_logLT()).
tls_surv_resid <- function(sp, target, p, model) {
  tdt_as_model_spec(model)$curve(sp$low, sp$up, sp$k * (target - sp$mid)) - p
}

# Root of a straight line through (ta, da) and (tb, db), interpolated on d =
# logLT - target when both ends are finite, else on the survival residual r
# (tls_surv_resid()); NA when neither is usable.
tls_interp_root <- function(da, db, ta, tb, ra = NULL, rb = NULL) {
  if (is.finite(da) && is.finite(db) && da != db) return(ta + da * (tb - ta) / (da - db))
  if (length(ra) && length(rb) && is.finite(ra) && is.finite(rb) && ra != rb)
    return(ta + ra * (tb - ta) / (ra - rb))
  NA_real_
}

# The single FALLING crossing in a row d = logLT - target (sign-safe: +-Inf
# allowed, NA undefined): the step where d goes from > 0 to <= 0 as temperature
# rises, i.e. where the LT shortens through t_ref. NA when there is none, or
# more than one: a rising crossing (the LT lengthening with temperature) only
# appears where a bent curve is extrapolated or near an asymptote, and two
# falling crossings leave CTmax ambiguous. `r` (optional) is the row's
# survival residual, used when d is infinite at an end of the step.
tls_falling_crossing <- function(d, grid, r = NULL) {
  n <- length(d)
  if (n < 2L) return(NA_real_)
  a <- d[-n]; b <- d[-1L]
  j <- which(!is.na(a) & !is.na(b) & a > 0 & b <= 0)
  if (length(j) != 1L) return(NA_real_)
  tls_interp_root(d[j], d[j + 1L], grid[j], grid[j + 1L], r[j], r[j + 1L])
}

# Which way to search beyond the grid for a draw with no crossing on it: +1
# (hotter) when the curve lies above the target across the grid, -1 (cooler)
# when it lies at or below it, 0 otherwise (undefined at that end, or a row
# with a rising crossing only). `dI` is the grid row of logLT - target
# (sign-safe).
tls_outward_side <- function(dI) {
  ok <- !is.na(dI); n <- length(dI)
  if (n && ok[n] && all(dI[ok] > 0))  return(1L)
  if (n && ok[1] && all(dI[ok] <= 0)) return(-1L)
  0L
}

# Per-draw crossing of `target`, searched first on the requested grid itself
# (columns `inner` of M, exactly as tls_invert_logLT()), then -- only for a draw
# with no crossing there whose curve lies wholly on one side of the target
# (tls_outward_side()) -- outward along the extension on the side the curve
# points to (above: hotter, since logLT falls with temperature; below: cooler),
# taking the first crossing beyond the grid: the one nearest the data. A bent
# curve extrapolated far from the data can turn back and cross the target
# again, so a plain inversion over the whole extended grid could return that
# spurious second crossing. The search stops (NA) at an undefined cell. `M` is
# sign-safe logLT (tls_clamp_logLT()); `R` its survival residuals, or NULL for
# a relative threshold. `grid` is sorted increasing (tls_extend_grid()) and
# `inner` the contiguous columns of the requested grid.
tls_invert_outward <- function(M, target, grid, inner = seq_along(grid), R = NULL) {
  Tc <- tls_invert_logLT(M[, inner, drop = FALSE], target, grid[inner],
                         R = if (is.null(R)) NULL else R[, inner, drop = FALSE])
  lo <- min(inner); hi <- max(inner); nc <- ncol(M)
  if (lo == 1L && hi == nc) return(Tc)
  for (i in which(!is.finite(Tc))) {
    side <- tls_outward_side(M[i, inner] - target)
    path <- if (side > 0L && hi < nc) c(hi, (hi + 1L):nc)
            else if (side < 0L && lo > 1L) c(lo, (lo - 1L):1L) else next
    d <- M[i, path] - target
    k <- if (side > 0L) which(is.na(d) | d <= 0)[1] else which(is.na(d) | d > 0)[1]
    if (is.na(k) || is.na(d[k])) next
    Tc[i] <- tls_interp_root(d[k - 1L], d[k], grid[path[k - 1L]], grid[path[k]],
                             if (is.null(R)) NULL else R[i, path[k - 1L]],
                             if (is.null(R)) NULL else R[i, path[k]])
  }
  Tc
}

# Second pass for bent LT curves, all groups together. Each element of `groups`
# is one group whose crossings were first inverted on its requested grid:
#   newrow    one-row data frame (moderators + filled columns);
#   tg        its sorted, unique centred grid;
#   sp        the tls_eval_subpars() list (low, up, k, mid, logLT) on the grid
#             columns, rows = the first pass's draws;
#   crossings named list of per-draw crossings (original scale; NA = none);
#   targets   matching named list of log10 reference durations (model units).
# Stage A refines the grid (tls_refine_grid(): points at most 0.05 C apart, or
# a twentieth of an interval wider than 1 C), evaluates the new points for
# every draw and inverts again on the refined grid, so a crossing is no longer
# bounded by the coarse grid's interpolation error. Stage B evaluates, for the
# draws still without a crossing whose curve points off the grid
# (tls_outward_side()), the refined grid's extension (tls_extend_grid(...,
# width)) and searches outward (tls_invert_outward()). posterior_linpred() has
# a large fixed cost per call, so in each stage the groups share engine calls;
# a group's draws are split into chunks when needed, so that no call
# evaluates more than `max_cells` draws x rows, and each chunk is inverted as
# soon as it is evaluated. Only the refined values of the draws Stage B needs
# are kept between the stages. `draw_ids` are the first pass's posterior draw
# ids (NULL = all, in order): row i of `sp` is posterior draw draw_ids[i].
# Returns `groups` with `tg` (refined) and the crossings updated, and `sp`
# dropped.
tls_extend_crossings <- function(fit, groups, temp, width, Tbar, draw_ids, bounds,
                                 nlpars = c("lowraw", "upraw", "logk", "mid"),
                                 re_formula = NA, mode = "absolute", p = NA,
                                 model = NULL, max_cells = 3e6) {
  pp   <- if (identical(mode, "absolute")) p else NA_real_
  nms  <- c(low = "low", up = "up", k = "k", mid = "mid", logLT = "logLT")
  ids_of <- function(rows) if (is.null(draw_ids)) rows else draw_ids[rows]
  resid  <- function(sp, tgt) if (is.na(pp)) NULL else tls_surv_resid(sp, tgt, pp, model)
  # TRUE for the rows of a clamped logLT matrix L whose crossing v is missing
  # and whose curve stays on one side of the target (the crossing may lie off
  # the grid, so Stage B searches outward for it).
  points_off <- function(L, v, target) {
    out <- logical(length(v))
    i <- which(!is.finite(v))
    out[i[vapply(i, function(r) tls_outward_side(L[r, ] - target) != 0L, logical(1))]] <- TRUE
    out
  }
  # Evaluate, for each group gi in `todo`, the temperatures plan[[gi]]$new at
  # the draws rows_of(gi) (first-pass row indices) and hand each block to
  # consume(gi, rows, sub), where `sub` is a tls_eval_subpars() list of
  # [length(rows) x length(new)] matrices. Groups share engine calls; a
  # group's draws are split into chunks so that no call exceeds `max_cells`.
  eval_batched <- function(plan, todo, rows_of, consume) {
    units <- list()
    for (gi in todo) {
      r <- rows_of(gi); nn <- length(plan[[gi]]$new)
      if (!length(r) || !nn) next
      size <- max(1L, floor(max_cells / nn))
      for (ch in split(r, ceiling(seq_along(r) / size)))
        units[[length(units) + 1L]] <- list(gi = gi, rows = ch)
    }
    cells <- function(us) length(unique(unlist(lapply(us, `[[`, "rows")))) *
      sum(vapply(us, function(u) length(plan[[u$gi]]$new), integer(1)))
    batches <- list(); cur <- list()
    for (u in units) {
      cand <- c(cur, list(u))
      if (length(cur) && cells(cand) > max_cells) {
        batches[[length(batches) + 1L]] <- cur; cur <- list(u)
      } else cur <- cand
    }
    if (length(cur)) batches[[length(batches) + 1L]] <- cur
    for (b in batches) {
      rows_u <- sort(unique(unlist(lapply(b, `[[`, "rows"))))
      nd <- do.call(rbind, lapply(b, function(u) {
        x <- groups[[u$gi]]$newrow[rep(1L, length(plan[[u$gi]]$new)), , drop = FALSE]
        x[[temp]] <- plan[[u$gi]]$new
        x
      }))
      rownames(nd) <- NULL
      spe <- tls_eval_subpars(fit, nd, bounds, nlpars = nlpars, re_formula = re_formula,
                              draw_ids = ids_of(rows_u), mode = mode,
                              p = if (is.na(pp)) 0.5 else pp, warn = FALSE)
      col0 <- 0L
      for (u in b) {
        nn   <- length(plan[[u$gi]]$new)
        cols <- col0 + seq_len(nn); col0 <- col0 + nn
        ru   <- match(u$rows, rows_u)
        consume(u$gi, u$rows, lapply(nms, function(nm) spe[[nm]][ru, cols, drop = FALSE]))
      }
    }
    invisible(NULL)
  }
  # `rows`: the first-pass rows whose values `sp` holds (all of them, until
  # Stage A keeps only the rows Stage B can need).
  for (gi in seq_along(groups)) groups[[gi]]$rows <- seq_len(nrow(groups[[gi]]$sp$logLT))

  # --- Stage A: refine the grid and invert again, every draw ---------------
  planA <- lapply(groups, function(g) {
    tf <- tls_refine_grid(g$tg)
    keep <- match(g$tg, tf)
    if (length(tf) == length(g$tg)) return(NULL)
    list(grid = tf, keep = keep, new = tf[-keep])
  })
  todoA <- which(!vapply(planA, is.null, logical(1)))
  if (length(todoA)) {
    kept <- vector("list", length(groups))
    consumeA <- function(gi, rows, sub) {
      pl <- planA[[gi]]; g <- groups[[gi]]
      spf <- lapply(nms, function(nm) {
        m <- matrix(NA_real_, length(rows), length(pl$grid))
        m[, pl$keep] <- g$sp[[nm]][rows, , drop = FALSE]
        m[, -pl$keep] <- sub[[nm]]
        m
      })
      L <- tls_clamp_logLT(spf, pp)
      off <- logical(length(rows))
      for (nm in names(g$crossings)) {
        v <- tls_invert_logLT(L, g$targets[[nm]], Tbar + pl$grid, R = resid(spf, g$targets[[nm]]))
        groups[[gi]]$crossings[[nm]][rows] <<- v
        off <- off | points_off(L, v, g$targets[[nm]])
      }
      if (any(off))
        kept[[gi]][[length(kept[[gi]]) + 1L]] <<-
          list(rows = rows[off], sp = lapply(spf, function(m) m[off, , drop = FALSE]))
    }
    eval_batched(planA, todoA, function(gi) groups[[gi]]$rows, consumeA)
    for (gi in todoA) {
      kk <- kept[[gi]]; ng <- length(planA[[gi]]$grid)
      groups[[gi]]$tg   <- planA[[gi]]$grid
      groups[[gi]]$rows <- as.integer(unlist(lapply(kk, `[[`, "rows")))
      groups[[gi]]$sp   <- lapply(nms, function(nm)
        if (length(kk)) do.call(rbind, lapply(kk, function(x) x$sp[[nm]]))
        else matrix(NA_real_, 0L, ng))
    }
    rm(kept)
  }

  # --- Stage B: search outward for the draws still without a crossing -------
  planB <- lapply(groups, function(g) {
    if (!length(g$rows)) return(NULL)
    L <- tls_clamp_logLT(g$sp, pp)
    off <- logical(length(g$rows))
    for (nm in names(g$crossings))
      off <- off | points_off(L, g$crossings[[nm]][g$rows], g$targets[[nm]])
    if (!any(off)) return(NULL)
    # a one-point grid has no spacing of its own: extend it at the finest step
    grid  <- tls_extend_grid(g$tg, width = width, max_per_side = 200L,
                             refine = if (length(g$tg) > 1L) 4L else 200L)
    inner <- match(g$tg, grid)
    list(local = which(off), need = g$rows[off], grid = grid, inner = inner,
         new = grid[-inner])
  })
  todoB <- which(!vapply(planB, is.null, logical(1)))
  if (length(todoB)) {
    consumeB <- function(gi, rows, sub) {
      pl <- planB[[gi]]; g <- groups[[gi]]
      loc <- pl$local[match(rows, pl$need)]
      spf <- lapply(nms, function(nm) {
        m <- matrix(NA_real_, length(rows), length(pl$grid))
        m[, pl$inner]  <- g$sp[[nm]][loc, , drop = FALSE]
        m[, -pl$inner] <- sub[[nm]]
        m
      })
      Lf <- tls_clamp_logLT(spf, pp)
      for (nm in names(g$crossings)) {
        v <- g$crossings[[nm]][rows]
        f <- which(!is.finite(v))
        if (!length(f)) next
        v[f] <- tls_invert_outward(
          Lf[f, , drop = FALSE], g$targets[[nm]], Tbar + pl$grid, pl$inner,
          R = resid(lapply(spf, function(m) m[f, , drop = FALSE]), g$targets[[nm]]))
        groups[[gi]]$crossings[[nm]][rows] <<- v
      }
    }
    eval_batched(planB, todoB, function(gi) planB[[gi]]$need, consumeB)
  }
  for (gi in seq_along(groups)) groups[[gi]]$sp <- NULL
  groups
}

# TRUE when every finite draw of a group's logLT is a straight line in
# temperature (tls_linear_draws()), so the closed-form crossing is exact; FALSE
# when any draw bends (or no draw is finite), so the crossing is inverted.
tls_group_is_linear <- function(M, temps) {
  lin <- tls_linear_draws(M, temps)
  any(!is.na(lin)) && all(lin, na.rm = TRUE)
}

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
# the disjoint-bounds 4PL). `M` is the evaluated logLT matrix
# `[ndraws x length(temp_grid)]`, `temp_grid` increasing; it may be sign-safe
# (tls_clamp_logLT(): +-Inf where the target lies outside the asymptotes),
# with `R` the matching survival residuals (tls_surv_resid()) to place a
# crossing next to such a cell. A row that is not monotone or not finite takes
# the single falling crossing (tls_falling_crossing()); it used to interpolate
# temperature against y-sorted values, which mixes the branches of a U- or
# cap-shaped row and can return a temperature that is not a crossing at all.
# Mirrors the retired tdt_ctmax_from_pars() for monotone rows.
tls_invert_logLT <- function(M, target, temp_grid, R = NULL) {
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
  for (i in which(bad))
    Tc[i] <- tls_falling_crossing(M[i, ] - target, temp_grid,
                                  if (is.null(R)) NULL else R[i, ])
  Tc
}

# Per-group z / CTmax / 1-hour-CTmax-anchor from logLT already evaluated on the
# group's grid. This is the SINGLE derivation contract used by tls() and the
# kept derive_* readers, so results are identical by construction: z is always the
# central-difference local z pooled over `z_grid` (tls_local_z; exact for a
# linear relative midpoint, accurate for a bent absolute curve), and the CTmax /
# 1-hour anchor use the exact closed form when the logLT curve is a straight
# line in temperature or a true per-draw inversion of the bent LT curve
# (tls_invert_logLT) otherwise. RNG-free -- the T_crit rate draws stay with the
# caller so each keeps its own reproducibility contract.
#
# logLT_plus / logLT_minus : [ndraws x length(z_grid)] logLT at z_grid +/- h.
# logLT_ctmax              : [ndraws x length(ctmax_grid)] logLT at ctmax_grid
#                            (`ctmax_grid` on the ORIGINAL temperature scale, so
#                            the returned CTmax is too; sorted increasing).
# mid0                     : [ndraws] logLT at temp_c = 0 (closed-form anchor).
# target / target_1hr      : log10 reference exposures on the model time scale.
# Tbar                     : temperature-centring constant (temp_mean).
# use_closed               : NA (default) decides from the curve itself
#                            (tls_group_is_linear() over every evaluated
#                            point); TRUE / FALSE force the closed form / the
#                            inversion. The parameterisation alone cannot
#                            decide: a direct fit's relative midpoint is
#                            linear, while a midpoint written with a
#                            non-linear temperature term is not.
# want_ctmax / want_ct1    : compute the CTmax / 1-hour anchor (skipped, left
#                            NULL, when a caller does not need them -- e.g.
#                            tls(params = "z")).
# sp_ctmax, p, model       : for an absolute threshold p, the tls_eval_subpars()
#                            list on the ctmax columns and the fit's model, so
#                            a bent curve is inverted sign-safely next to cells
#                            where p lies outside the asymptotes
#                            (tls_clamp_logLT(), tls_surv_resid()). p = NA
#                            (relative) inverts logLT_ctmax as is.
# Returns list(z = <tls_local_z output>, ctmax = <per-draw vector or NULL>,
#              ct1 = <per-draw vector or NULL>, closed = <logical or NA when
#              neither CTmax nor the anchor was asked for>).
tls_zct_draws <- function(logLT_plus, logLT_minus, h, z_grid,
                          logLT_ctmax, ctmax_grid, mid0,
                          target, target_1hr, Tbar, use_closed = NA,
                          z_local = FALSE, want_ctmax = TRUE, want_ct1 = FALSE,
                          probs = c(0.025, 0.5, 0.975), sp_ctmax = NULL,
                          p = NA, model = NULL) {
  z_obj <- tls_local_z(logLT_plus, logLT_minus, h, z_grid,
                       probs = probs, local = isTRUE(z_local))
  ctmax <- NULL; ct1 <- NULL; closed <- NA
  if (isTRUE(want_ctmax) || isTRUE(want_ct1)) {
    Tb <- Tbar %||% 0
    closed <- if (is.na(use_closed)) {
      tls_group_is_linear(cbind(mid0, logLT_minus, logLT_plus, logLT_ctmax),
                          c(0, z_grid - h, z_grid + h, ctmax_grid - Tb))
    } else isTRUE(use_closed)
    if (closed) {
      # Straight logLT(T) (e.g. any linear relative midpoint, midpoint or
      # direct): the t_ref crossing is closed-form at any temperature.
      # `slope` is the mean central-difference slope, which equals
      # b_mid_temp_c for a linear midpoint.
      slope <- rowMeans((logLT_plus - logLT_minus) / (2 * h))
      if (isTRUE(want_ctmax)) ctmax <- Tb + (target     - mid0) / slope
      if (isTRUE(want_ct1))   ct1   <- Tb + (target_1hr - mid0) / slope
    } else {
      # Bent curve (an absolute threshold whose asymptotes or steepness change
      # with temperature): invert the fitted logLT surface per draw at the
      # reference exposure on the grid. Draws that do not cross there are left
      # NA for the caller's outward search (tls_extend_crossings()).
      abs_p <- !is.null(sp_ctmax) && !is.na(p)
      L  <- if (abs_p) tls_clamp_logLT(sp_ctmax, p) else logLT_ctmax
      Rf <- function(tgt) if (abs_p) tls_surv_resid(sp_ctmax, tgt, p, model) else NULL
      if (isTRUE(want_ctmax)) ctmax <- tls_invert_logLT(L, target,     ctmax_grid, R = Rf(target))
      if (isTRUE(want_ct1))   ct1   <- tls_invert_logLT(L, target_1hr, ctmax_grid, R = Rf(target_1hr))
    }
  }
  list(z = z_obj, ctmax = ctmax, ct1 = ct1, closed = closed)
}
