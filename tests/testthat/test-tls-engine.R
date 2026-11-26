# Fast, brms-free tests for the engine plumbing in R/tls_engine.R that does not
# call posterior_linpred: the grid builder and the by-resolution. (tls_local_z /
# tls_invert_logLT are covered in test-tdt-z-local-gate.R and
# test-tdt-ctmax-vectorised.R; tls_eval_subpars / tls_draw_ids are exercised
# gated through every migrated reader on the cached fixtures.)

test_that("tls_build_grid: ungrouped grid is temp_grid long, fills model columns, .grp = 'all'", {
  mdata <- data.frame(temp_c = c(-2, 0, 2), logd = 1, life_stage = factor("a"))
  nd <- tls_build_grid(mdata, by = NULL, temp = "temp_c", temp_grid = c(-1, 0, 1))
  expect_equal(nd$temp_c, c(-1, 0, 1))
  expect_true(all(nd$.grp == "all"))
  expect_true("life_stage" %in% names(nd))           # missing model column filled from row 1
  expect_true(all(nd$life_stage == "a"))
})

test_that("tls_build_grid: grouped grid crosses moderator levels x temp_grid, tagged in .grp", {
  mdata <- data.frame(temp_c = c(-2, 0, 2),
                      life_stage = factor(c("young", "old", "young")))
  nd <- tls_build_grid(mdata, by = "life_stage", temp = "temp_c", temp_grid = c(0, 1))
  expect_setequal(unique(nd$.grp), c("young", "old"))
  expect_equal(nrow(nd), 2L * 2L)                     # 2 levels x 2 temps
  expect_equal(sum(nd$life_stage == "young"), 2L)
  expect_setequal(unique(nd$temp_c), c(0, 1))
})

test_that("tdt_resolve_by: explicit wins; NULL uses meta$group_vars; single-condition -> NULL", {
  wf_g <- structure(list(meta = list(group_vars = "sex")), class = "bayes_tls")
  wf_s <- structure(list(meta = list(group_vars = character(0))), class = "bayes_tls")
  expect_equal(tdt_resolve_by(wf_g, by = NULL), "sex")       # auto from the fit
  expect_equal(tdt_resolve_by(wf_g, by = "clone"), "clone")  # explicit overrides
  expect_null(tdt_resolve_by(wf_s, by = NULL))               # single-condition
})

test_that("tdt_resolve_by safety net: grouped-looking fit with no recorded moderators errors", {
  set.seed(1); np <- 30
  df <- data.frame(b_CTmaxdev_Intercept = rnorm(np), b_CTmaxdev_grpB = rnorm(np),
                   b_logz_Intercept = rnorm(np))
  wf <- structure(list(fit = posterior::as_draws_df(df), meta = list()),  # no group_vars
                  class = "bayes_tls")
  expect_error(tdt_resolve_by(wf, by = NULL), "by =")        # don't silently use the reference level
  expect_equal(tdt_resolve_by(wf, by = "grp"), "grp")        # explicit by still works
})

# --- Draw subsamples, bounds and crossings (2026-11-26) ----------------------

test_that("tls_draw_ids() validates ndraws before touching the fit", {
  expect_null(tls_draw_ids(NULL, NULL))                # NULL = all draws
  for (bad in list(0, -5, NA, NA_real_, 2.5, Inf, "10", c(10, 20), TRUE))
    expect_error(tls_draw_ids(NULL, bad), "single whole number >= 1")
  expect_silent(tdt_check_ndraws(100L))
  expect_silent(tdt_check_ndraws(100))
})

test_that("tdt_workflow_bounds() uses meta$bounds, else rebuilds them from lower/upper", {
  b <- compute_4pl_bounds(0, 1)
  expect_identical(tdt_workflow_bounds(list(meta = list(bounds = "as fitted"))), "as fitted")
  expect_equal(tdt_workflow_bounds(list(meta = list())), b)
  expect_equal(tdt_workflow_bounds(list(meta = list(lower = 0.1, upper = 0.9))),
               compute_4pl_bounds(0.1, 0.9))
})

test_that("tls_eval_subpars() refuses incomplete bounds before calling brms", {
  expect_error(tls_eval_subpars(NULL, data.frame(), list(low_min = 0)),
               "compute_4pl_bounds")
  expect_error(tls_eval_subpars(NULL, data.frame(), NULL), "compute_4pl_bounds")
})

test_that("tls_abs_outside_message() counts draws and gives the share of evaluations", {
  bad <- matrix(FALSE, 1000, 6); bad[1, 1] <- bad[2, 4] <- bad[2, 5] <- TRUE
  low <- matrix(0.03, 1000, 6)
  up  <- matrix(0.97, 1000, 6); up[1, 1] <- 0.69; up[2, 4:5] <- 0.6928
  m <- tls_abs_outside_message(bad, low, up, 0.7)
  expect_match(m, "outside the fitted curve's asymptotes for 2 of 1000 draws", fixed = TRUE)
  expect_match(m, "(3 of 6000 evaluations, 0.05%)", fixed = TRUE)
  expect_match(m, "at or below 0.7 (survival starts below the target) in 3", fixed = TRUE)
  expect_match(m, "(survival never falls to it) in 0", fixed = TRUE)
  m2 <- tls_abs_outside_message(matrix(TRUE, 2, 2), matrix(0.6, 2, 2),
                                matrix(0.9, 2, 2), 0.5)
  expect_match(m2, "for 2 of 2 draws", fixed = TRUE)
  expect_match(m2, "(4 of 4 evaluations, 100%)", fixed = TRUE)
  expect_match(m2, "never falls to it) in 4", fixed = TRUE)
  m3 <- tls_abs_outside_message(matrix(TRUE, 1, 1), matrix(NA_real_, 1, 1),
                                matrix(NA_real_, 1, 1), 0.5)
  expect_match(m3, "1 have undefined asymptotes", fixed = TRUE)
})

test_that("tls_linear_draws() separates straight from bent curves per draw", {
  tt <- c(3, -3, 0, -1.5, 1.5, 0)                       # unsorted, with a repeat
  M  <- rbind(2 - 0.3 * tt,                             # straight
              2 - 0.3 * tt + 0.01 * tt^2,               # bent
              2 - 0.3 * tt + 1e-12 * tt^2,              # straight to rounding
              rep(1, length(tt)),                       # flat (straight)
              c(NA, 2 - 0.3 * tt[-1]))                  # non-finite
  expect_identical(tls_linear_draws(M, tt), c(TRUE, FALSE, TRUE, TRUE, NA))
  expect_true(all(is.na(tls_linear_draws(M, rep(0, length(tt))))))   # no span
  expect_true(tls_group_is_linear(M[c(1, 3, 5), ], tt))  # NA draws ignored
  expect_false(tls_group_is_linear(M[1:2, ], tt))        # any bent draw: invert
  expect_false(tls_group_is_linear(M[5, , drop = FALSE], tt))   # nothing finite
})

test_that("tls_extend_grid() keeps the grid and extends it finely on both sides", {
  tg <- seq(-3, 3, length.out = 11)
  g  <- tls_extend_grid(tg, width = 6)
  expect_equal(g[match(tg, g)], tg)                     # grid points kept exactly
  expect_false(anyNA(match(tg, g)))
  expect_equal(range(g), c(-9, 9))
  expect_equal(unique(round(diff(g[g <= -3]), 10)), 0.15)   # a quarter of the spacing
  expect_equal(unique(round(diff(g[g >= 3]), 10)), 0.15)
  expect_true(all(diff(g) > 0))
  expect_equal(range(tls_extend_grid(tg, width = 0)), c(-9, 9))  # at least the grid width
  expect_lte(length(tls_extend_grid(tg, width = 1000)), 11 + 2 * 80)
  expect_equal(tls_extend_grid(0), seq(-1, 1, by = 0.25))   # single point: 1 degree a side
})

test_that("tls_invert_outward() finds the first crossing beyond the grid", {
  grid  <- seq(-9, 9, by = 0.5)
  inner <- which(grid >= -3 & grid <= 3)
  lin   <- function(a, b) a + b * grid                  # straight logLT(T)
  # inside the grid: exactly tls_invert_logLT() on the grid
  M <- rbind(lin(0, -0.3))
  expect_identical(tls_invert_outward(M, 0.2, grid, inner),
                   tls_invert_logLT(M[, inner, drop = FALSE], 0.2, grid[inner]))
  # wholly above the target on the grid -> hotter; wholly below -> cooler
  M <- rbind(lin(2, -0.3), lin(-2, -0.3))
  expect_equal(tls_invert_outward(M, 0, grid, inner), c(2 / 0.3, -2 / 0.3))
  expect_true(all(is.na(tls_invert_logLT(M[, inner], 0, grid[inner]))))
  # a curve that turns back before reaching the target stays NA
  tb <- ifelse(grid <= 5, 2 - 0.2 * grid, 1 + 0.5 * (grid - 5))
  expect_true(is.na(tls_invert_outward(rbind(tb), 0, grid, inner)))
  # a curve that crosses, then turns back and crosses again: the first crossing
  y2 <- ifelse(grid <= 6, 2 - 0.4 * grid, -0.4 + 1.2 * (grid - 6))
  got <- tls_invert_outward(rbind(y2), 0, grid, inner)
  expect_equal(got, 5)
  # a non-finite value before the crossing stops the search
  y3 <- lin(2, -0.3); y3[grid == 4] <- NA
  expect_true(is.na(tls_invert_outward(rbind(y3), 0, grid, inner)))
  # every column inner: identical to tls_invert_logLT()
  M <- rbind(lin(2, -0.3), y2)
  expect_identical(tls_invert_outward(M, 0, grid), tls_invert_logLT(M, 0, grid))
})


test_that("tls_zct_draws() uses the closed form for a straight curve, beyond the grid too", {
  set.seed(3)
  np <- 50; zg <- seq(-3, 3, length.out = 7); h <- 1e-3
  a  <- rnorm(np, 1, 0.1); b <- -runif(np, 0.1, 0.3)  # logLT = a + b * temp_c
  f  <- function(tt) outer(a, rep(1, length(tt))) + outer(b, tt)
  target <- 2.5                                         # crossing far below the grid
  zc <- tls_zct_draws(f(zg + h), f(zg - h), h, zg, f(zg), zg + 33, a,
                      target = target, target_1hr = 0, Tbar = 33,
                      want_ctmax = TRUE, want_ct1 = TRUE)
  expect_true(zc$closed)
  expect_equal(zc$ctmax, 33 + (target - a) / b)
  expect_true(all(zc$ctmax < 33 - 3))                   # outside the grid, still exact
  # A bent curve is inverted on the grid only (NA beyond it, for the outward pass)
  g  <- function(tt) f(tt) + 0.02 * outer(rep(1, np), tt^2)
  zb <- tls_zct_draws(g(zg + h), g(zg - h), h, zg, g(zg), zg + 33, g(0),
                      target = target, target_1hr = 0, Tbar = 33)
  expect_false(zb$closed)
  expect_true(all(is.na(zb$ctmax)))
  # forcing the inversion on the straight curve reproduces the grid behaviour
  zf <- tls_zct_draws(f(zg + h), f(zg - h), h, zg, f(zg), zg + 33, a,
                      target = target, target_1hr = 0, Tbar = 33, use_closed = FALSE)
  expect_true(all(is.na(zf$ctmax)))
})

test_that("tls_check_temp_grid() catches degrees C grids that miss the centred data", {
  # snow gum: 30-56 degrees C, temp_mean 42.64 -> centred -12.64 to 13.36
  sg <- c(-12.64, 13.36)
  expect_warning(tls_check_temp_grid(seq(30, 38, 2), sg, 42.64), "entirely outside")
  expect_warning(tls_check_temp_grid(seq(30, 38, 2), sg, 42.64), "temp_grid - 42.64")
  # aphid cold branch: centred -3.417 to 3.583, temp_mean -7.58
  expect_warning(tls_check_temp_grid(c(-10, -7, -4), c(-3.417, 3.583), -7.58),
                 "entirely outside")
  # centred grids overlapping the data stay silent, even if they extend beyond it
  expect_silent(tls_check_temp_grid(c(-3, 9), c(-3, 3), 33))
  expect_silent(tls_check_temp_grid(seq(-20, 20, 5), sg, 42.64))
  expect_warning(tls_check_temp_grid(c(-12, 0), c(-3, 0, 3), 33), "reaches far outside")
})

# --- Sign-safe crossings near the asymptotes (second review, 2026-11-26) ------

test_that("tdt_finite_grid() drops non-finite temperatures with a warning", {
  expect_null(tdt_finite_grid(NULL))
  expect_identical(tdt_finite_grid(c(1, 2)), c(1, 2))
  expect_warning(g <- tdt_finite_grid(c(1, NA, Inf, 2, -Inf)), "Dropping 3 non-finite")
  expect_identical(g, c(1, 2))
  expect_error(tdt_finite_grid(c(NA, Inf)), "no finite temperatures")
})

test_that("tls_clamp_logLT() and tls_surv_resid() give the sign of logLT - target everywhere", {
  set.seed(5)
  n <- 200
  low <- matrix(stats::runif(n, 0, 0.45), n, 1); up <- matrix(stats::runif(n, 0.55, 1), n, 1)
  k   <- matrix(stats::rexp(n) + 0.5, n, 1);     mid <- matrix(stats::rnorm(n), n, 1)
  for (model in c("log-logistic", "log-normal", "weibull1", "weibull2")) {
    for (p in c(0.1, 0.5, 0.9)) {
      sp <- list(low = low, up = up, k = k, mid = mid,
                 logLT = mid + tdt_abs_offset(model, low, up, k, p))
      L  <- tls_clamp_logLT(sp, p)
      expect_true(all(!is.na(L)), info = paste(model, p))
      expect_true(all(L[up <= p] == -Inf) && all(L[low >= p] == Inf), info = paste(model, p))
      for (tgt in c(-1, 0, 1)) {
        r <- tls_surv_resid(sp, tgt, p, model)
        expect_true(all(is.finite(r)))
        expect_equal(sign(r), sign(L - tgt), info = paste(model, p, tgt))
      }
    }
  }
  sp_rel <- list(logLT = matrix(c(1, NA), 1), low = matrix(0, 1, 2), up = matrix(1, 1, 2))
  expect_identical(tls_clamp_logLT(sp_rel, NA), sp_rel$logLT)   # relative: unchanged
})

test_that("tls_falling_crossing() takes the single crossing where the LT shortens", {
  # U shape (the old y-sorted fallback returned 32.93, which is not a crossing)
  u_row <- c(0.0593, 0.0220, -0.0106, -0.0373, -0.0563, -0.0651, -0.0603, -0.0372,
             0.0115, 0.0974, 0.2408)
  g <- seq(30, 36, by = 0.6)
  expect_equal(tls_falling_crossing(u_row, g), 30.6 + 0.6 * 0.0220 / 0.0326)
  expect_equal(tls_invert_logLT(rbind(u_row), 0, g), 30.6 + 0.6 * 0.0220 / 0.0326)
  # cap shape (a real snow gum draw at p = 0.9: the old fallback gave 30.9)
  cap <- c(1.80, 2.25, 2.35, 2.21, 1.91, 1.50, 1.00, 0.45, -0.14, -0.76, -1.40)
  gs  <- seq(30, 56, length.out = 11)
  expect_equal(tls_invert_logLT(rbind(cap), log10(60), gs),
               gs[5] + (1.91 - log10(60)) * (gs[6] - gs[5]) / (1.91 - 1.50))
  # no falling crossing, a rising one only, or two falling ones: NA
  expect_true(is.na(tls_falling_crossing(c(1, 2, 3), 1:3)))
  expect_true(is.na(tls_falling_crossing(c(-1, 0.5, 2), 1:3)))
  expect_true(is.na(tls_falling_crossing(c(1, -1, 1, -1), 1:4)))
  # a step next to an infinite LT is placed on the survival residual
  expect_equal(tls_falling_crossing(c(0.5, -Inf), c(10, 11), r = c(0.02, -0.06)), 10.25)
  expect_equal(tls_falling_crossing(c(Inf, -0.5), c(10, 11), r = c(0.06, -0.02)), 10.75)
  expect_true(is.na(tls_falling_crossing(c(0.5, -Inf), c(10, 11))))       # no residual
  expect_equal(tls_invert_logLT(rbind(c(1, 0.5, -Inf, -Inf)), 0, 1:4,
                                R = rbind(c(0.3, 0.1, -0.1, -0.3))), 2.5)
})

test_that("tls_outward_side() points the search off the grid only when the row allows it", {
  expect_identical(tls_outward_side(c(3, 2, 1)), 1L)              # above: hotter
  expect_identical(tls_outward_side(c(-1, -2, -Inf)), -1L)        # at/below: cooler
  expect_identical(tls_outward_side(c(0, -1, -2)), -1L)
  expect_identical(tls_outward_side(c(-1, 1, 2)), 0L)             # rising only
  expect_identical(tls_outward_side(c(3, 2, NA)), 0L)             # undefined hot end
  expect_identical(tls_outward_side(c(NA, -1, -2)), 0L)           # undefined cool end
  expect_identical(tls_outward_side(c(Inf, 2, NA, 1)), 1L)        # interior NA ignored
})

test_that("tls_invert_outward() crosses next to an undefined cell and from partly undefined rows", {
  grid  <- seq(-9, 9, by = 0.5)
  inner <- which(grid >= -3 & grid <= 3)
  # above the target on the grid, then the upper asymptote falls to p:
  # logLT -> -Inf; the crossing lies in that last step, placed on the residual
  y <- 2 - 0.3 * grid; y[grid > 5] <- -Inf               # y(5) = 0.5 > 0.2
  r <- 0.21 - 0.04 * grid                                # residual crosses 0 at 5.25
  expect_equal(tls_invert_outward(rbind(y), 0.2, grid, inner, R = rbind(r)), 5.25)
  expect_true(is.na(tls_invert_outward(rbind(y), 0.2, grid, inner)))   # no residual
  # below the target on the grid, the hot end undefined (-Inf): search cooler
  y2 <- -0.5 - 0.3 * grid; y2[grid > 2] <- -Inf
  expect_equal(tls_invert_outward(rbind(y2), 0, grid, inner), -0.5 / 0.3)
  # a genuinely undefined cell (NA) on the path still stops the search
  y3 <- 2 - 0.3 * grid; y3[grid == 4] <- NA
  expect_true(is.na(tls_invert_outward(rbind(y3), 0, grid, inner)))
  # the first crossing along the path is taken even after the curve turns
  # away from the target (it never did on the fits checked)
  y4 <- ifelse(grid <= 4, 2 - 0.3 * grid, ifelse(grid <= 6, 0.8 + 0.2 * (grid - 4),
                                                 1.2 - (grid - 6)))
  expect_equal(tls_invert_outward(rbind(y4), 0, grid, inner), 7.2)
})

test_that("readers that pass ndraws on to brms validate it first", {
  wf0 <- structure(list(fit = structure(list(), class = "brmsfit")), class = "bayes_tls")
  expect_error(posterior_linpred_tdt(wf0, data.frame(), ndraws = 2.5), "whole number")
  tr <- data.frame(time = c(0, 1), temp = c(30, 30))
  expect_error(predict_heat_injury(tr, wf0, ndraws = 0), "whole number")
})

test_that("tls_refine_grid() keeps the grid points and fills each step to <= 0.05 C", {
  tg <- seq(-3, 3, length.out = 11)                      # 0.6 C steps
  g  <- tls_refine_grid(tg)
  expect_false(anyNA(match(tg, g)))                     # original points kept exactly
  expect_equal(length(g), 121L)
  expect_lt(max(diff(g)), 0.05 + 1e-12)
  expect_true(all(diff(g) > 0))
  wide <- tls_refine_grid(seq(-13, 13, length.out = 11)) # 2.6 C steps: at most 20 a step
  expect_equal(length(wide), 10L * 20L + 1L)
  expect_identical(tls_refine_grid(seq(0, 1, by = 0.05)), seq(0, 1, by = 0.05))  # fine already
  expect_identical(tls_refine_grid(c(2, 0, 2)), tls_refine_grid(c(0, 2)))   # sorted, unique
  expect_equal(length(tls_refine_grid(c(0, 2))), 21L)   # a 2 C step: capped at 20 (0.1 C)
  expect_identical(tls_refine_grid(1), 1)
})

# --- tls_extend_crossings() on an analytic engine (no brms) -------------------
# tls_eval_subpars() is replaced by an exact 4PL whose draws are rows of `par`,
# so every crossing has an exact reference (uniroot on the survival curve).

mock_engine <- function(par, mid_fun, seen = NULL) {
  force(par); force(mid_fun)
  function(fit, newdata, bounds, nlpars = NULL, re_formula = NA, draw_ids = NULL,
           mode = "absolute", p = 0.5, warn = TRUE) {
    i  <- if (is.null(draw_ids)) seq_len(nrow(par)) else draw_ids
    tt <- newdata$temp_c
    if (!is.null(seen)) seen$calls[[length(seen$calls) + 1L]] <- list(t = tt, n = length(i))
    low <- matrix(0.05, length(i), length(tt))
    up  <- pmin(par$u0[i] + outer(par$us[i], tt), 0.999)
    k   <- matrix(par$k[i], length(i), length(tt))
    mid <- mid_fun(par[i, , drop = FALSE], tt)
    logLT <- if (identical(mode, "absolute")) mid + tdt_abs_offset("log-logistic", low, up, k, p) else mid
    list(low = low, up = up, k = k, mid = mid, logLT = logLT)
  }
}

# First pass as tls() builds it: invert each draw on the coarse grid `tg`.
first_pass <- function(engine, tg, target, p, rows = NULL) {
  sp <- engine(NULL, data.frame(temp_c = tg), NULL, draw_ids = rows,
               mode = if (is.na(p)) "relative" else "absolute", p = if (is.na(p)) 0.5 else p)
  R  <- if (is.na(p)) NULL else tls_surv_resid(sp, target, p, "log-logistic")
  list(newrow = data.frame(temp_c = 0), tg = tg, sp = sp, targets = list(v = target),
       crossings = list(v = tls_invert_logLT(tls_clamp_logLT(sp, p), target, tg, R = R)))
}

# Exact crossing per draw: the single falling crossing of S(t_ref) - p (or of
# logLT - target) inside [lo, hi], else the first one outward up to `reach`.
exact_crossings <- function(f, n, lo, hi, reach) {
  fine <- seq(lo - reach, hi + reach, by = 0.001)
  ins  <- fine >= lo - 1e-9 & fine <= hi + 1e-9
  vapply(seq_len(n), function(d) {
    y <- f(d, fine); yi <- y[ins]; fi <- fine[ins]
    j <- which(yi[-length(yi)] > 0 & yi[-1] <= 0)
    root <- function(a, b) stats::uniroot(function(t) f(d, t), c(a, b), tol = 1e-10)$root
    if (length(j) == 1L) return(root(fi[j], fi[j + 1L]))
    if (length(j) > 1L) return(NA_real_)
    if (all(yi > 0)) { o <- which(fine > hi); k <- which(y[o] <= 0)[1]
      return(if (is.na(k)) NA_real_ else root(fine[o][k] - 0.001, fine[o][k])) }
    if (all(yi <= 0)) { o <- rev(which(fine < lo)); k <- which(y[o] > 0)[1]
      return(if (is.na(k)) NA_real_ else root(fine[o][k], fine[o][k] + 0.001)) }
    NA_real_
  }, numeric(1))
}

test_that("tls_extend_crossings() refines on the grid and steps outward on the refined grid", {
  set.seed(42)
  par <- data.frame(u0 = runif(200, 0.84, 0.91), us = -runif(200, 0.015, 0.03),
                    m0 = runif(200, 1.2, 1.6), ms = -runif(200, 0.15, 0.25), k = runif(200, 4, 8))
  midf <- function(q, tt) q$m0 + outer(q$ms, tt)
  seen <- new.env(); seen$calls <- list()
  eng <- mock_engine(par, midf, seen)
  local_mocked_bindings(tls_eval_subpars = eng)
  p <- 0.8; target <- 0; tg <- seq(-3, 3, length.out = 11)     # coarse: 0.6 C steps
  g <- first_pass(eng, tg, target, p)
  seen$calls <- list()
  out <- tls_extend_crossings(NULL, list(a = g), "temp_c", width = 6, Tbar = 0,
                              draw_ids = NULL, bounds = NULL, mode = "absolute", p = p,
                              model = "log-logistic")$a
  surv <- function(d, t) {
    up <- pmin(par$u0[d] + par$us[d] * t, 0.999)
    0.05 + (up - 0.05) / (1 + exp(par$k[d] * (target - (par$m0[d] + par$ms[d] * t)))) - p
  }
  ref <- exact_crossings(surv, 200, -3, 3, 6)
  v   <- out$crossings$v
  expect_identical(is.finite(v), is.finite(ref))
  on  <- is.finite(ref) & abs(ref) <= 3; off <- is.finite(ref) & abs(ref) > 3
  expect_gt(sum(on), 30); expect_gt(sum(off), 30)
  expect_lt(max(abs(v[on] - ref[on])), 0.05)      # a refined step (coarse grid: about 0.15)
  expect_lt(max(abs(v[off] - ref[off])), 0.03)    # a refined outward step (coarse: about 0.05)
  expect_equal(out$tg, tls_refine_grid(tg))
  beyond <- sort(unique(unlist(lapply(seen$calls, `[[`, "t"))))
  beyond <- beyond[beyond > max(tg)]
  expect_gt(length(beyond), 0)
  expect_lte(max(diff(c(max(tg), beyond))), 6 / 200 + 1e-9)   # the refined grid's extension
  expect_null(out$sp)
})

test_that("tls_extend_crossings() handles a relative bent midpoint, mixed groups, chunks and draw ids", {
  set.seed(7)
  par <- data.frame(u0 = 0.95, us = 0, m0 = runif(150, 1.0, 1.4), ms = -runif(150, 0.15, 0.25),
                    mq = -runif(150, 0.002, 0.01), k = 6)
  midf <- function(q, tt) q$m0 + outer(q$ms, tt) + outer(q$mq, tt^2)
  seen <- new.env(); seen$calls <- list()
  eng <- mock_engine(par, midf, seen)
  local_mocked_bindings(tls_eval_subpars = eng)
  target <- 0.5
  ref <- exact_crossings(function(d, t) par$m0[d] + par$ms[d] * t + par$mq[d] * t^2 - target,
                         150, -3, 3, 6)
  coarse <- seq(-3, 3, length.out = 11)
  fine   <- seq(-3, 3, by = 0.05)
  run <- function(groups, ...) tls_extend_crossings(NULL, groups, "temp_c", width = 6, Tbar = 0,
                                                    draw_ids = NULL, bounds = NULL,
                                                    mode = "relative", p = NA, ...)
  both <- run(list(c = first_pass(eng, coarse, target, NA), f = first_pass(eng, fine, target, NA)))
  ok <- is.finite(ref)
  expect_gt(sum(ok), 100)
  for (g in c("c", "f")) {
    expect_identical(is.finite(both[[g]]$crossings$v), ok)
    expect_lt(max(abs(both[[g]]$crossings$v[ok] - ref[ok])), 1e-3)
  }
  expect_equal(both$f$tg, fine)                    # already fine: left as it is
  expect_equal(both$c$tg, tls_refine_grid(coarse))
  alone <- run(list(c = first_pass(eng, coarse, target, NA)))
  expect_identical(both$c$crossings, alone$c$crossings)
  # draws split into chunks give the same crossings, and a subsample keeps its draw ids
  g1 <- first_pass(eng, coarse, target, NA)
  seen$calls <- list()
  chunked <- run(list(c = g1), max_cells = 500)
  expect_identical(chunked$c$crossings, alone$c$crossings)
  cells <- vapply(seen$calls, function(x) x$n * length(x$t), numeric(1))
  expect_gt(length(cells), 1)
  expect_lte(max(cells), 500)                      # no engine call exceeds max_cells
  ids <- sort(sample.int(150, 40))
  sub <- tls_extend_crossings(NULL, list(c = first_pass(eng, coarse, target, NA, rows = ids)),
                              "temp_c", width = 6, Tbar = 0, draw_ids = ids, bounds = NULL,
                              mode = "relative", p = NA, max_cells = 500)
  expect_equal(sub$c$crossings$v, alone$c$crossings$v[ids])
})
