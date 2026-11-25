# Exact, fit-free tests of predict_heat_injury(shape = "varying"). The real
# function body runs on known synthetic curves; only the accessors that read a
# brms fit are mocked. They check both methods against the references in
# helper-heat-injury.R (every model) and against mechanistic truths, their
# documented properties (for the per-individual method: no order or start
# effect without repair, monotone in added time, heat injury of rank p reaching
# 100% when survival falls to p; for both: the start value, permanent deaths on
# exposure, equality with the constant shape when the shape is flat), and
# branches the cached fixtures cannot reach (a lower asymptote that changes with
# temperature, every model's far tail, dropped draws, repeated time stamps,
# grouped fits, a model fitted in minutes). They run without RUN_BRMS_TESTS.

# Per-draw curve pieces as functions of temperature, each returning a
# draws x temperatures matrix.
lin_T <- function(at_tbar, slope, tbar = 33) {
  function(temps) {
    outer(at_tbar, rep(1, length(temps))) +
      outer(rep_len(slope, length(at_tbar)), temps - tbar)
  }
}
exp_T <- function(at_tbar, rate, tbar = 33) {
  function(temps) outer(at_tbar, exp(rate * (temps - tbar)))
}

# A stand-in workflow whose curves are `spec`: one list per group with `low`,
# `up`, `k` (functions of temperature), `m0`, `m1` (per draw; mid is linear in
# temperature, mid(T) = m0 + m1 (T - tbar), as the constant shape assumes) and
# an optional `bar` overriding the centring-temperature parameters that the
# constant shape and the repair conversion read. `drop_local` lists, per group,
# draws whose local curve counts as invalid at a trace temperature.
local_synthetic_wf <- function(spec, model = "log-logistic", tbar = 33,
                               drop_local = list(), unit = "hours",
                               env = parent.frame()) {
  grouped <- !is.null(names(spec))
  labels  <- if (grouped) names(spec) else "all"
  pars <- do.call(rbind, lapply(seq_along(spec), function(g) {
    s <- spec[[g]]
    b <- list(low = s$low(tbar)[, 1], up = s$up(tbar)[, 1], k = s$k(tbar)[, 1])
    if (!is.null(s$bar)) b <- utils::modifyList(b, s$bar)
    out <- data.frame(.draw = seq_along(s$m0), low = b$low, up = b$up, k = b$k,
                      mid_int = s$m0, mid_temp = s$m1)
    if (grouped) out$grp <- labels[g]
    out
  }))
  local_curves <- function(temps) {
    out <- lapply(seq_along(spec), function(g) {
      s    <- spec[[g]]
      keep <- setdiff(seq_along(s$m0), drop_local[[labels[g]]])
      mid  <- lin_T(s$m0, s$m1, tbar)(temps)
      list(temps = temps, .draw = keep, grp = labels[g],
           low = s$low(temps)[keep, , drop = FALSE],
           up  = s$up(temps)[keep, , drop = FALSE],
           k   = s$k(temps)[keep, , drop = FALSE],
           mid = mid[keep, , drop = FALSE],
           by_row = if (grouped) data.frame(grp = labels[g]) else NULL)
    })
    stats::setNames(out, labels)
  }
  testthat::local_mocked_bindings(
    has_fit         = function(workflow) TRUE,
    tdt_model       = function(x) model,
    tdt_resolve_by  = function(workflow, by) if (grouped) "grp" else NULL,
    hi_pars         = function(workflow, by = NULL) pars,
    hi_local_curves = function(workflow, by = NULL, temps_unique) local_curves(temps_unique),
    .env = env)
  list(meta = list(temp_mean = tbar, duration_unit = unit))
}

# Per-trace-point curve values of one draw.
draw_P <- function(spec_g, d, temps, tbar = 33) {
  data.frame(low = spec_g$low(temps)[d, ], up = spec_g$up(temps)[d, ],
             k = spec_g$k(temps)[d, ],
             mid = lin_T(spec_g$m0, spec_g$m1, tbar)(temps)[d, ])
}

repair_pars_test <- function(r_ref) {
  list(TA = 14065, TAL = 50000, TAH = 100000, TL = 15 + 273.15,
       TH = 35 + 273.15, TREF = 25 + 273.15, r_ref = r_ref)
}
rep_rate_of <- function(temp, rp) {
  repair_rate_schoolfield(temp, TA = rp$TA, TAL = rp$TAL, TAH = rp$TAH,
                          TL = rp$TL, TH = rp$TH, TREF = rp$TREF, r_ref = rp$r_ref)
}
# Rows ordered by draw, then time: the draw order of separate calls differs.
by_draw <- function(d) d[order(d$.draw, d$time), ]
end_rows <- function(o) {
  d <- o$draws[o$draws$time == max(o$draws$time), ]
  d[order(d$.draw), ]
}

# Three draws whose low, up and k all change with temperature.
spec_moving <- list(list(
  low = lin_T(c(0.04, 0.05, 0.06), 0.01),
  up  = lin_T(c(0.91, 0.92, 0.93), -0.015),
  k   = exp_T(c(5.5, 6, 6.5), 0.08),
  m0  = c(1.38, 1.40, 1.42), m1 = c(-0.17, -0.18, -0.19)))

test_that("varying shape follows the per-individual rule when low, up and k change with temperature", {
  # The trace starts at 30 C, below T_c, so survival starts at the 32 C top.
  temps <- c(30, 32, 34, 36)
  tr <- data.frame(time = 0:72)
  tr$temp <- temps[tr$time %% 4 + 1]
  rp <- repair_pars_test(0.05)
  rr <- rep_rate_of(tr$temp, rp)
  s  <- spec_moving[[1]]
  for (m in c("log-logistic", "log-normal", "weibull1", "weibull2")) {
    wf <- local_synthetic_wf(spec_moving, model = m)
    for (rep_on in c(FALSE, TRUE)) {
      for (tg in list("relative", 0.5)) {
        for (irr in c(TRUE, FALSE)) {
          o <- predict_heat_injury(tr, wf, target_surv = tg, T_c = 30.5, repair = rep_on,
                                   repair_pars = if (rep_on) rp,
                                   irreversible_mortality = irr, shape = "varying",
                                   save_draws = TRUE)
          for (d in 1:3) {
            conv <- 1; p <- NULL
            if (is.numeric(tg)) {   # r_ref in doses of the target, at T_bar, by the model
              p <- tg
              conv <- 10^tdt_model_spec(m)$abs_corr(s$low(33)[d, 1], s$up(33)[d, 1],
                                                    s$k(33)[d, 1], p)
            }
            want <- ref_individual(tr$time, tr$temp, draw_P(s, d, tr$temp), model = m,
                                   p_target = p, conv = conv, T_c = 30.5,
                                   rep_rate = if (rep_on) rr, irreversible = irr)
            got  <- o$draws[o$draws$.draw == d, ]
            lab  <- paste(m, "repair", rep_on, "target", tg, "irreversible", irr, "draw", d)
            expect_equal(got$dose, want$dose, tolerance = 1e-10, info = lab)
            expect_equal(got$survival, want$survival, tolerance = 1e-10, info = lab)
          }
        }
      }
    }
  }
})

test_that("the surviving rank is exact when low rises with temperature (kinks in injury)", {
  # A lower asymptote that rises with temperature gives each individual's
  # injury kinks with an infinite slope where its rank crosses a bin's low, and
  # the root often sits just above one. A solver that stopped on a short Newton
  # step there returned the previous survival (errors near 1e-3 on this trace).
  spec_kink <- list(list(low = lin_T(c(0.0867, 0.09), 0.0169),
                         up  = lin_T(c(0.812, 0.82), -0.0065),
                         k   = exp_T(c(13.5, 11), 0.0138),
                         m0  = c(1.17, 1.2), m1 = c(-0.231, -0.231)))
  tm <- 0:336
  tr <- data.frame(time = tm, temp = round(33 + 1.95 * sin(2 * pi * tm / 24 + 2), 2))
  for (m in c("log-normal", "weibull1")) {
    wf <- local_synthetic_wf(spec_kink, model = m)
    o  <- predict_heat_injury(tr, wf, irreversible_mortality = FALSE, shape = "varying",
                              save_draws = TRUE)$draws
    for (d in 1:2) {
      want <- ref_individual(tr$time, tr$temp, draw_P(spec_kink[[1]], d, tr$temp),
                             model = m, irreversible = FALSE)
      expect_equal(o$survival[o$.draw == d], want$survival, tolerance = 1e-10,
                   info = paste(m, d))
    }
  }
})

test_that("every model: exact at a constant temperature and after a change of temperature", {
  # A shallow hot curve leaves survival moderate while the dose is far past
  # the steep cool curve's midpoint (deep in that curve's tail).
  spec_tail <- list(list(low = lin_T(rep(0.05, 3), 0), up = lin_T(rep(0.9, 3), 0),
                         k = exp_T(c(15.1, 15.4, 15.7), -0.68),
                         m0 = rep(1.4, 3), m1 = rep(-0.18, 3)))
  tr  <- data.frame(time = 0:90, temp = c(rep(36, 60), rep(31, 31)))
  hot <- data.frame(time = 0:60, temp = 36)
  for (m in c("log-logistic", "log-normal", "weibull1", "weibull2")) {
    wf <- local_synthetic_wf(spec_tail, model = m)
    o  <- predict_heat_injury(tr, wf, shape = "varying", save_draws = TRUE)
    for (d in 1:3) {
      want <- ref_individual(tr$time, tr$temp, draw_P(spec_tail[[1]], d, tr$temp), model = m)
      got  <- o$draws[o$draws$.draw == d, ]
      expect_true(all(is.finite(got$survival)), info = m)
      expect_equal(got$survival, want$survival, tolerance = 1e-10, info = paste(m, d))
    }
    oc <- predict_heat_injury(hot, wf, shape = "varying", save_draws = TRUE)$draws
    P  <- draw_P(spec_tail[[1]], 1, hot$temp)
    curve_1 <- tdt_model_spec(m)$curve(P$low, P$up,
                                      P$k * log10(pmax(hot$time / 10^P$mid, 1e-12)))
    expect_equal(oc$survival[oc$.draw == 1], curve_1, tolerance = 1e-10, info = m)
  }
})

test_that("individuals with their own TDT lines and fixed ranks: survival is exact", {
  # Individual Z has log10 LT(T) = a0 + b0 dT + Z (sa + sb dT), Z ~ N(0, 1), and
  # dies when sum(dt / LT) reaches 1. Ranks never cross, so the per-individual
  # rule must reproduce the population's survival under any trace; the truth
  # is solved for the critical Z independently. At constant T the curves are
  # the log-normal model with low 0, up 1, mid a0 + b0 dT, k 1 / (sa + sb dT).
  a0 <- 1.4; b0 <- -0.18; sa <- 0.2; sb <- 0.04
  spec <- list(list(low = lin_T(c(0, 0), 0), up = lin_T(c(1, 1), 0),
                    k = function(temps) outer(c(1, 1), 1 / (sa + sb * (temps - 33))),
                    m0 = c(a0, a0), m1 = c(b0, b0)))
  wf <- local_synthetic_wf(spec, model = "log-normal")
  tm <- seq(0, 72, by = 0.25)
  tr <- data.frame(time = tm, temp = round(33 + 3 * sin(2 * pi * tm / 24), 4))
  truth <- numeric(nrow(tr)); truth[1] <- 1
  dT <- tr$temp - 33; dt <- diff(tr$time)
  for (j in 2:nrow(tr)) {
    f <- function(z) log(sum(dt[1:(j - 1)] * 10^(-(a0 + b0 * dT[1:(j - 1)]) -
                                                  z * (sa + sb * dT[1:(j - 1)]))))
    truth[j] <- stats::pnorm(stats::uniroot(f, c(-40, 40), tol = 1e-14)$root,
                             lower.tail = FALSE)
  }
  o <- predict_heat_injury(tr, wf, shape = "varying", save_draws = TRUE)$draws
  expect_equal(o$survival[o$.draw == 1], truth, tolerance = 1e-9)
  expect_lt(tail(truth, 1), 0.5)               # the trace kills most of them
})

test_that("without repair, neither the order of temperatures nor the start of the trace matters", {
  wf <- local_synthetic_wf(spec_moving)
  hot_first  <- data.frame(time = (0:144) / 6, temp = c(rep(36, 72), rep(31, 73)))
  cool_first <- data.frame(time = (0:144) / 6, temp = c(rep(31, 72), rep(36, 73)))
  end_s <- function(tr) end_rows(predict_heat_injury(tr, wf, shape = "varying",
                                                     save_draws = TRUE))$survival
  expect_equal(end_s(hot_first), end_s(cool_first), tolerance = 1e-12)
  # whole daily cycles, logged from the coolest or the hottest hour
  cyc <- function(sg) {
    tm <- seq(0, 96, by = 0.25)
    data.frame(time = tm, temp = round(33 + sg * 3 * cos(2 * pi * tm / 24), 4))
  }
  expect_equal(end_s(cyc(-1)), end_s(cyc(1)), tolerance = 1e-12)
})

test_that("extra heat never raises survival", {
  wf <- local_synthetic_wf(spec_moving)
  week     <- data.frame(time = 0:168, temp = 31)
  week_hot <- data.frame(time = 0:171, temp = c(36, 36, 36, rep(31, 169)))
  afternoons <- data.frame(time = 0:168,
                           temp = ifelse((0:168) %% 24 >= 12 & (0:168) %% 24 < 15, 36, 31))
  s_w  <- end_rows(predict_heat_injury(week, wf, shape = "varying", save_draws = TRUE))
  s_wh <- end_rows(predict_heat_injury(week_hot, wf, shape = "varying", save_draws = TRUE))
  s_af <- end_rows(predict_heat_injury(afternoons, wf, shape = "varying", save_draws = TRUE))
  expect_true(all(s_wh$survival < s_w$survival))
  expect_true(all(s_af$survival < s_w$survival))
})

test_that("heat injury of an absolute target reaches 100% exactly when survival falls to it", {
  wf <- local_synthetic_wf(spec_moving)
  tm <- seq(0, 96, by = 0.25)
  tr <- data.frame(time = tm, temp = round(33 + 3 * sin(2 * pi * tm / 24), 4))
  o  <- predict_heat_injury(tr, wf, target_surv = 0.6, shape = "varying",
                            save_draws = TRUE)$draws
  away <- abs(o$survival - 0.6) > 1e-9          # skip a step landing exactly on 0.6
  expect_true(any(o$hi >= 100) && any(o$hi < 100))
  expect_identical((o$hi >= 100)[away], (o$survival <= 0.6)[away])
})

test_that("varying shape starts at up(T_1) and never falls below the floor", {
  wf <- local_synthetic_wf(spec_moving)
  tr <- data.frame(time = 0:300, temp = c(36, rep(31, 300)))
  o  <- predict_heat_injury(tr, wf, shape = "varying", save_draws = TRUE)
  s  <- spec_moving[[1]]
  first <- o$draws[o$draws$time == 0, ]
  expect_equal(first$survival[order(first$.draw)], as.vector(s$up(36)), tolerance = 1e-9)
  expect_true(all(end_rows(o)$survival >= as.vector(s$low(31)) - 1e-12))
})

# Two draws whose top rises with temperature: a cool start has the lowest top.
spec_up <- list(list(low = lin_T(rep(0.05, 2), 0), up = lin_T(c(0.80, 0.82), 0.03),
                     k = lin_T(c(6, 5.5), 0), m0 = c(1.4, 1.38), m1 = rep(-0.18, 2)))

test_that("a start that accrues no injury sets nothing: survival starts where injury does", {
  wf <- local_synthetic_wf(spec_up)            # up(20) is 0.41, up(36) 0.89
  s  <- spec_up[[1]]
  hot        <- data.frame(time = 0:12, temp = 36)
  cool_start <- data.frame(time = 0:22, temp = c(rep(20, 10), rep(36, 13)))
  dup_start  <- data.frame(time = c(0, 0:12), temp = c(20, rep(36, 13)))
  for (meth in c("individual", "cumulative")) {
    run <- function(tr, T_c = NULL) {
      by_draw(predict_heat_injury(tr, wf, T_c = T_c, shape = "varying",
                                  varying_method = meth, save_draws = TRUE)$draws)
    }
    h <- run(hot)
    # 10 h at 20 C, below T_c: survival starts at the 36 C top and the hot hours
    # then give exactly what they give alone
    cs <- run(cool_start, T_c = 30)
    expect_equal(cs$survival[cs$time == 0], as.vector(s$up(36)), tolerance = 1e-9, info = meth)
    expect_equal(cs$survival[cs$time >= 10], h$survival, tolerance = 1e-12, info = meth)
    # a repeated first time stamp at 20 C is no exposure either (no T_c)
    ds <- run(dup_start)
    first <- !duplicated(ds$.draw)
    expect_equal(ds$survival[first], h$survival[h$time == 0], tolerance = 1e-12, info = meth)
    expect_equal(ds$survival[!first], h$survival, tolerance = 1e-12, info = meth)
    # a trace that accrues no injury at all keeps the constant shape's value
    none <- data.frame(time = 0:6, temp = c(20, 25, 36, 36, 22, 20, 20))
    nv <- run(none, T_c = 40)
    nc <- by_draw(predict_heat_injury(none, wf, T_c = 40, save_draws = TRUE)$draws)
    expect_equal(nv$survival, nc$survival, tolerance = 1e-12, info = meth)
    expect_true(all(nv$hi == 0), info = meth)
  }
})

test_that("from a cool start below T_c, order still does not matter and added time never raises survival", {
  wf <- local_synthetic_wf(spec_up)
  end_s <- function(tr, irr = TRUE) {
    end_rows(predict_heat_injury(tr, wf, T_c = 30, irreversible_mortality = irr,
                                 shape = "varying", save_draws = TRUE))$survival
  }
  a <- data.frame(time = 0:30, temp = c(rep(20, 6), rep(31, 12), rep(36, 13)))
  b <- data.frame(time = 0:30, temp = c(rep(20, 6), rep(36, 12), rep(31, 13)))
  expect_equal(end_s(a), end_s(b), tolerance = 1e-12)
  expect_equal(end_s(a, irr = FALSE), end_s(a), tolerance = 1e-12)   # no repair
  # mild enough that survival stays well above the 20 C top (0.41): a start
  # value read at 20 C would pin the base trace below the one with more heat
  base <- data.frame(time = 0:16, temp = c(rep(20, 10), rep(33, 7)))
  plus <- data.frame(time = 0:19, temp = c(36, 36, 36, rep(20, 10), rep(33, 7)))
  expect_true(all(end_s(base) > 0.6))
  expect_true(all(end_s(plus) < end_s(base)))
})

test_that("deaths on exposure to a lower top are permanent, also when mortality is reversible", {
  # up falls with temperature: 5 h at 31 C, 5 h at 38 C, then a cool spell
  # whose repair clears all injury. Survival rebounds to the 38 C top (0.80),
  # not to the 31 C top it started from (0.94).
  spec_dn <- list(list(low = lin_T(rep(0.05, 2), 0), up = lin_T(c(0.9, 0.91), -0.02),
                       k = lin_T(c(6, 5.5), 0), m0 = c(1.4, 1.38), m1 = rep(-0.18, 2)))
  wf <- local_synthetic_wf(spec_dn)
  s  <- spec_dn[[1]]
  tr <- data.frame(time = 0:60, temp = c(rep(31, 5), rep(38, 5), rep(20, 51)))
  for (meth in c("individual", "cumulative")) {
    d <- by_draw(predict_heat_injury(tr, wf, T_c = 25, repair = TRUE,
                                     repair_pars = repair_pars_test(0.2),
                                     irreversible_mortality = FALSE, shape = "varying",
                                     varying_method = meth, save_draws = TRUE)$draws)
    for (k in 1:2) {
      dk <- d[d$.draw == k, ]
      top38 <- s$up(38)[k, 1]
      expect_true(all(dk$survival[dk$time >= 6] <= top38 + 1e-12), info = meth)
      expect_lt(min(dk$survival), 0.5 * top38)                         # injury killed many
      expect_equal(dk$hi[nrow(dk)], 0)                                 # all repaired
      expect_equal(tail(dk$survival, 1), top38, tolerance = 1e-9, info = meth)
    }
  }
})

test_that("with reversible mortality, survival rebounds from deep in the tail (flat shape = constant)", {
  # A steep curve driven 36 doses deep: the at-risk share is ~1e-20, below the
  # rounding of survival near low. Repair then removes the dose, and survival
  # must rebound exactly as the constant shape's g(dose) does.
  flat <- list(list(low = lin_T(rep(0.05, 2), 0), up = lin_T(rep(0.9, 2), 0),
                    k = lin_T(c(30, 28), 0), m0 = c(1.4, 1.38), m1 = c(-0.18, -0.18)))
  wf <- local_synthetic_wf(flat)
  tr <- data.frame(time = 0:250, temp = c(rep(40, 50), rep(20, 201)))
  rp <- repair_pars_test(1)
  args <- list(tr, wf, T_c = 30, repair = TRUE, repair_pars = rp,
               repair_scales_with_survival = FALSE, irreversible_mortality = FALSE,
               save_draws = TRUE)
  hc <- by_draw(do.call(predict_heat_injury, c(args, shape = "constant"))$draws)
  hv <- by_draw(do.call(predict_heat_injury, c(args, shape = "varying"))$draws)
  # reached the floor (to the solver's precision there, about 1e-13 in survival)
  expect_lt(min(hv$survival), 0.05 + 1e-12)
  expect_gt(max(end_rows(list(draws = hv))$survival), 0.5)   # and rebounded
  expect_equal(hv$survival, hc$survival, tolerance = 1e-10)
  expect_equal(hv$hi, hc$hi, tolerance = 1e-10)
})

test_that("a repeated time stamp adds nothing, even where the target is already passed", {
  spec <- list(list(low = lin_T(rep(0.05, 2), 0), up = lin_T(c(0.95, 0.96), -0.03),
                    k = lin_T(rep(6, 2), 0), m0 = c(1.4, 1.41), m1 = rep(-0.18, 2)))
  wf <- local_synthetic_wf(spec)               # up(37) is about 0.83, below 0.85
  once  <- data.frame(time = c(0, 1, 2, 3, 4), temp = c(33, 37, 33, 33, 33))
  twice <- data.frame(time = c(0, 1, 1, 2, 3, 4), temp = c(33, 37, 37, 33, 33, 33))
  run <- function(tr) {
    expect_warning(o <- predict_heat_injury(tr, wf, target_surv = 0.85,
                                            shape = "varying", save_draws = TRUE),
                   "already passed")
    by_draw(o$draws)
  }
  a <- run(once); b <- run(twice)
  expect_false(anyNA(b$hi))
  expect_true(all(is.infinite(b$hi[b$time >= 2])))
  # the repeated row is the row before it, and every other row is unchanged
  dup <- which(duplicated(b[c(".draw", "time")]))
  expect_identical(b$survival[dup], b$survival[dup - 1])
  expect_identical(b$hi[dup], b$hi[dup - 1])
  expect_equal(b$survival[-dup], a$survival, tolerance = 1e-12)
  expect_equal(b$hi[-dup], a$hi, tolerance = 1e-12)
})

test_that("a passed target warns only where its rate is used, and T_c silences it", {
  spec <- list(list(low = lin_T(rep(0.05, 2), 0), up = lin_T(c(0.95, 0.96), -0.03),
                    k = lin_T(rep(6, 2), 0), m0 = c(1.4, 1.41), m1 = rep(-0.18, 2)))
  wf <- local_synthetic_wf(spec)
  last_only <- data.frame(time = 0:3, temp = c(33, 33, 33, 37))
  expect_no_warning(o <- predict_heat_injury(last_only, wf, target_surv = 0.85,
                                             shape = "varying", save_draws = TRUE))
  expect_true(all(is.finite(o$draws$hi)))
  hot_below_tc <- data.frame(time = 0:4, temp = c(33, 37, 37, 33, 33))
  expect_no_warning(o2 <- predict_heat_injury(hot_below_tc, wf, target_surv = 0.85,
                                              T_c = 37, shape = "varying",
                                              save_draws = TRUE))
  d1 <- o2$draws[o2$draws$.draw == 1, ]
  expect_equal(d1$hi[3], d1$hi[2])             # the 37 C interval adds nothing
  expect_true(all(is.finite(o2$draws$hi)))
})

test_that("a target below a temperature's low adds no heat injury there and leaves survival alone", {
  spec <- list(list(low = lin_T(c(0.2, 0.21), 0.05), up = lin_T(rep(0.95, 2), 0),
                    k = lin_T(rep(6, 2), 0), m0 = c(1.4, 1.41), m1 = rep(-0.18, 2)))
  wf <- local_synthetic_wf(spec)               # low(36) is 0.35, above 0.3
  tr <- data.frame(time = 0:6, temp = c(33, 36, 36, 33, 33, 36, 33))
  ab  <- by_draw(predict_heat_injury(tr, wf, target_surv = 0.3, shape = "varying",
                                     save_draws = TRUE)$draws)
  rel <- by_draw(predict_heat_injury(tr, wf, shape = "varying", save_draws = TRUE)$draws)
  expect_equal(ab$survival, rel$survival, tolerance = 1e-12)
  d1 <- ab[ab$.draw == 1, ]
  expect_equal(diff(d1$hi)[tr$temp[-nrow(tr)] == 36], c(0, 0, 0))
  expect_true(all(diff(d1$hi)[tr$temp[-nrow(tr)] == 33] > 0))
})

test_that("dropped draws are blanked, counted in meta$ndraws, and reported", {
  tr <- data.frame(time = 0:12, temp = rep(c(33, 35), length.out = 13))
  # (a) a local curve invalid at a trace temperature
  wf <- local_synthetic_wf(spec_moving, drop_local = list(all = 2))
  expect_warning(o <- predict_heat_injury(tr, wf, shape = "varying", save_draws = TRUE),
                 "1 of 3 posterior draws gave no prediction.*invalid local curve")
  expect_equal(o$meta$ndraws, 2L)
  expect_true(all(is.na(o$draws$survival[o$draws$.draw == 2])))
  expect_true(all(is.na(o$draws$hi[o$draws$.draw == 2])))
  # (b) repair with an absolute target outside a draw's curve at T_bar
  sp <- spec_moving
  sp[[1]]$bar <- list(up = c(0.40, 0.92, 0.93))
  wf <- local_synthetic_wf(sp)
  expect_warning(o <- predict_heat_injury(tr, wf, target_surv = 0.5, repair = TRUE,
                                          repair_pars = repair_pars_test(0.05),
                                          shape = "varying", save_draws = TRUE),
                 "could not be converted")
  expect_equal(o$meta$ndraws, 2L)
  expect_true(all(is.na(o$draws$hi[o$draws$.draw == 1])))
  # (c) every draw dropped: an NA summary and a warning, not an error
  wf <- local_synthetic_wf(spec_moving, drop_local = list(all = 1:3))
  expect_warning(o <- predict_heat_injury(tr, wf, shape = "varying"),
                 "3 of 3 posterior draws gave no prediction")
  expect_equal(o$meta$ndraws, 0L)
  expect_true(all(is.na(o$summary$surv_median)))
  expect_equal(nrow(o$summary), nrow(tr))
})

test_that("the constant shape leaves a dropped draw out of the heat-injury summary too", {
  sp <- spec_moving
  sp[[1]]$bar <- list(up = c(0.80, 0.92, 0.93))   # 0.85 unreachable at T_bar on draw 1
  wf <- local_synthetic_wf(sp)
  tr <- data.frame(time = 0:12, temp = 34)
  expect_warning(o <- predict_heat_injury(tr, wf, target_surv = 0.85,
                                          repair_scales_with_survival = FALSE,
                                          save_draws = TRUE),
                 "gave no prediction")
  expect_equal(o$meta$ndraws, 2L)
  expect_true(all(is.na(o$draws$hi[o$draws$.draw == 1])))
  kept <- o$draws[o$draws$.draw != 1 & o$draws$time == 12, ]
  expect_equal(tail(o$summary$hi_median, 1), stats::median(kept$hi))
})

test_that("grouped fits report meta$ndraws as the smallest group count, with each group's", {
  spec2 <- list(A = spec_moving[[1]], B = spec_moving[[1]])
  tr <- data.frame(time = 0:12, temp = rep(c(33, 35), length.out = 13))
  wf <- local_synthetic_wf(spec2, drop_local = list(B = 2))
  expect_warning(o <- predict_heat_injury(tr, wf, shape = "varying"), "gave no prediction")
  expect_equal(o$meta$ndraws, 2L)
  expect_equal(o$meta$ndraws_by_group, c(A = 3L, B = 2L))
  wf <- local_synthetic_wf(spec2, drop_local = list(B = 1:3))
  expect_warning(o <- predict_heat_injury(tr, wf, shape = "varying"),
                 "No draw gave a prediction for B")
  expect_equal(o$meta$ndraws_by_group, c(A = 3L, B = 0L))
  expect_true(all(is.na(o$summary$surv_median[o$summary$grp == "B"])))
  expect_true(all(is.finite(o$summary$surv_median[o$summary$grp == "A"])))
  # nothing dropped: one count, the same for every group
  wf <- local_synthetic_wf(spec2)
  expect_no_warning(o <- predict_heat_injury(tr, wf, shape = "varying"))
  expect_equal(o$meta$ndraws, 3L)
})

test_that("a grouped fit integrates each group with its own curves and repair conversion", {
  spec_a <- spec_moving[[1]]
  spec_b <- list(low = lin_T(c(0.08, 0.07, 0.09), 0.005), up = lin_T(c(0.85, 0.88, 0.9), -0.01),
                 k = exp_T(c(4, 4.5, 5), 0.05), m0 = c(1.2, 1.25, 1.3),
                 m1 = c(-0.15, -0.16, -0.14))
  tr <- data.frame(time = 0:48, temp = rep(c(31, 33, 35, 36), length.out = 49))
  args <- list(target_surv = 0.5, repair = TRUE, repair_pars = repair_pars_test(0.1),
               shape = "varying", save_draws = TRUE)
  for (meth in c("individual", "cumulative")) {
    wf_g <- local_synthetic_wf(list(A = spec_a, B = spec_b))
    g <- do.call(predict_heat_injury, c(list(tr, wf_g, varying_method = meth), args))$draws
    for (grp in c("A", "B")) {
      wf_1 <- local_synthetic_wf(list(if (grp == "A") spec_a else spec_b))
      one  <- by_draw(do.call(predict_heat_injury,
                              c(list(tr, wf_1, varying_method = meth), args))$draws)
      gg   <- by_draw(g[g$grp == grp, ])
      expect_equal(gg$survival, one$survival, tolerance = 1e-12, info = paste(meth, grp))
      expect_equal(gg$hi, one$hi, tolerance = 1e-12, info = paste(meth, grp))
    }
  }
})

test_that("a model fitted in minutes gives the same varying-shape prediction as in hours", {
  sp_min <- spec_moving
  sp_min[[1]]$m0 <- spec_moving[[1]]$m0 + log10(60)     # mid in log10 minutes
  tr <- data.frame(time = 0:48, temp = rep(c(31, 33, 35, 36), length.out = 49))
  args <- list(target_surv = 0.5, repair = TRUE, repair_pars = repair_pars_test(0.05),
               shape = "varying", save_draws = TRUE)
  for (meth in c("individual", "cumulative")) {
    wf_h <- local_synthetic_wf(spec_moving)
    h <- by_draw(do.call(predict_heat_injury, c(list(tr, wf_h, varying_method = meth), args))$draws)
    wf_m <- local_synthetic_wf(sp_min, unit = "minutes")
    m <- by_draw(do.call(predict_heat_injury, c(list(tr, wf_m, varying_method = meth), args))$draws)
    expect_equal(m$survival, h$survival, tolerance = 1e-10, info = meth)
    expect_equal(m$hi, h$hi, tolerance = 1e-10, info = meth)
  }
})

# ---- varying_method = "cumulative" (cumulative exposure) ------------------------

test_that("cumulative exposure follows its definition when low, up and k change with temperature", {
  temps <- c(30, 32, 34, 36)
  tr <- data.frame(time = 0:72)
  tr$temp <- temps[tr$time %% 4 + 1]
  rp <- repair_pars_test(0.05)
  rr <- rep_rate_of(tr$temp, rp)
  s  <- spec_moving[[1]]
  for (m in c("log-logistic", "log-normal", "weibull1", "weibull2")) {
    wf <- local_synthetic_wf(spec_moving, model = m)
    for (rep_on in c(FALSE, TRUE)) {
      for (irr in c(TRUE, FALSE)) {
        o <- predict_heat_injury(tr, wf, T_c = 30.5, repair = rep_on,
                                 repair_pars = if (rep_on) rp,
                                 irreversible_mortality = irr, shape = "varying",
                                 varying_method = "cumulative", save_draws = TRUE)
        for (d in 1:3) {
          want <- ref_cumulative(tr$time, tr$temp, draw_P(s, d, tr$temp), model = m,
                                 T_c = 30.5, rep_rate = if (rep_on) rr, irreversible = irr)
          lab <- paste(m, "repair", rep_on, "irreversible", irr, "draw", d)
          expect_equal(o$draws$survival[o$draws$.draw == d], want, tolerance = 1e-9,
                       info = lab)
        }
      }
    }
  }
  expect_identical(o$meta$varying_method, "cumulative")
})

test_that("cumulative exposure keeps its dose when one asymptote is the same at both temperatures", {
  # Close to an asymptote the gap to it underflows on the linear scale (the
  # Weibull type 2 curve near its top, the type 1 curve near its floor), so the
  # state must be carried in logs. Closed forms on alternating temperatures:
  # with equal tops, Weibull type 2 keeps log(up - S) = log(r) - log(2) e^(-u);
  # with equal floors, Weibull type 1 keeps log(S - low) = log(r) - log(2) e^u.
  tm <- seq(0, 96, by = 0.25)
  tr <- data.frame(time = tm, temp = ifelse(seq_along(tm) %% 2 == 1, 31, 32))
  closed <- function(P, type) {
    n <- nrow(P); r <- P$up - P$low
    surv <- function(u, j) {
      if (type == 2) P$up[j] - r[j] * 2^(-exp(-u)) else P$low[j] + r[j] * 2^(-exp(u))
    }
    S <- numeric(n); S[1] <- surv(P$k[1] * log10(1e-12), 1); u <- -Inf
    for (j in 2:n) {
      a <- j - 1
      if (a > 1 && r[a] != r[a - 1]) {          # carry the state to this step's curve
        x <- (if (type == 2) exp(-u) else exp(u)) + log(r[a] / r[a - 1]) / log(2)
        u <- if (type == 2) -log(x) else if (x > 0) log(x) else -Inf
      }
      u <- P$k[a] * log10(10^(u / P$k[a]) + diff(tm)[a] / 10^P$mid[a])
      S[j] <- min(S[a], surv(u, a))
    }
    S
  }
  # Weibull type 2: up 1 or 0.999 at both temperatures, low 0.04 / 0.05
  top_eq <- list(list(low = lin_T(c(0.06, 0.06), 0.01), up = lin_T(c(1, 0.999), 0),
                      k = lin_T(c(6, 6), 0), m0 = c(1.4, 1.4), m1 = c(-0.18, -0.18)))
  # Weibull type 1: low 0.05 at both, up 0.95 / 0.90; deep past the floor
  # (about 20 doses), where the excess above the floor underflows
  bot_eq <- list(list(low = lin_T(c(0.05, 0.05), 0), up = lin_T(c(0.85, 0.85), -0.05),
                      k = lin_T(c(6, 6), 0), m0 = c(0.4, 0.4), m1 = c(-0.18, -0.18)))
  for (cs in list(list(spec = top_eq, model = "weibull2", type = 2),
                  list(spec = bot_eq, model = "weibull1", type = 1))) {
    wf <- local_synthetic_wf(cs$spec, model = cs$model)
    o  <- by_draw(predict_heat_injury(tr, wf, shape = "varying", varying_method = "cumulative",
                                      save_draws = TRUE)$draws)
    for (d in 1:2) {
      want <- closed(draw_P(cs$spec[[1]], d, tr$temp), cs$type)
      expect_equal(o$survival[o$.draw == d], want, tolerance = 1e-10,
                   info = paste(cs$model, d))
    }
    expect_lt(max(o$survival[o$time == max(tm)]), 0.5)       # it does accrue injury
  }
})

test_that("cumulative exposure: a state below a temperature's floor stays on its curve, and repair still acts", {
  # low is 0.05 at 31 C and 0.55 at 36 C: after 8 h at 31 C survival (about
  # 0.32) is below the 36 C floor, so the 36 C hour adds no injury to it, and
  # injury continues at 31 C afterwards (not frozen).
  spec_fl <- list(list(low = function(temps) outer(c(1, 1), ifelse(temps >= 36, 0.55, 0.05)),
                       up = lin_T(c(0.95, 0.95), 0), k = lin_T(c(6, 5.5), 0),
                       m0 = c(0.4, 0.42), m1 = c(-0.18, -0.18)))
  wf <- local_synthetic_wf(spec_fl)
  tr <- data.frame(time = 0:15, temp = c(rep(31, 8), 36, rep(31, 7)))
  rp <- repair_pars_test(0.02); rr <- rep_rate_of(tr$temp, rp)
  for (rep_on in c(FALSE, TRUE)) {
    for (irr in c(TRUE, FALSE)) {
      o <- by_draw(predict_heat_injury(tr, wf, repair = rep_on, repair_pars = if (rep_on) rp,
                                       repair_scales_with_survival = FALSE,
                                       irreversible_mortality = irr, shape = "varying",
                                       varying_method = "cumulative", save_draws = TRUE)$draws)
      for (d in 1:2) {
        want <- ref_cumulative(tr$time, tr$temp, draw_P(spec_fl[[1]], d, tr$temp),
                               rep_rate = if (rep_on) rr, scale = FALSE, irreversible = irr)
        got  <- o$survival[o$.draw == d]
        lab  <- paste("repair", rep_on, "irreversible", irr, "draw", d)
        expect_equal(got, want, tolerance = 1e-10, info = lab)
        expect_lt(got[9], 0.55)                                   # below the 36 C floor
        expect_lt(tail(got, 1), got[10] - 0.05)                   # and still falling after
      }
    }
  }
  # Weibull type 1, low the same everywhere: 110 h at 36 C drives survival far
  # past the floor; a cool spell of strong repair must still bring it back to
  # the 36 C top with reversible mortality.
  spec_w1 <- list(list(low = lin_T(c(0.05, 0.05), 0),
                       up = function(temps) outer(c(1, 1), ifelse(temps >= 36, 0.85, 0.95)),
                       k = lin_T(c(6, 6), 0), m0 = c(1.4, 1.4), m1 = c(-0.18, -0.18)))
  wf <- local_synthetic_wf(spec_w1, model = "weibull1")
  tr <- data.frame(time = 0:260, temp = c(rep(36, 110), rep(20, 151)))
  rp <- list(TA = 14065, TAL = 50000, TAH = 100000, TL = 15 + 273.15,
             TH = 22 + 273.15, TREF = 20 + 273.15, r_ref = 1)
  for (meth in c("individual", "cumulative")) {
    d1 <- by_draw(predict_heat_injury(tr, wf, repair = TRUE, repair_pars = rp,
                                      repair_scales_with_survival = FALSE,
                                      irreversible_mortality = FALSE, shape = "varying",
                                      varying_method = meth, save_draws = TRUE)$draws)
    expect_lt(min(d1$survival), 0.051)
    expect_equal(d1$hi[nrow(d1)], 0)
    expect_equal(d1$survival[d1$time == 260], c(0.85, 0.85), tolerance = 1e-9, info = meth)
  }
})

test_that("cumulative exposure: injury partly repaired below T_c follows the definition", {
  tr <- data.frame(time = 0:40, temp = c(rep(34, 12), rep(25, 22), rep(34, 7)))
  rp <- repair_pars_test(0.005); rr <- rep_rate_of(tr$temp, rp)
  s  <- spec_moving[[1]]
  wf <- local_synthetic_wf(spec_moving)
  for (scl in c(TRUE, FALSE)) {
    for (irr in c(TRUE, FALSE)) {
      args <- list(tr, wf, T_c = 26, repair = TRUE, repair_pars = rp,
                   repair_scales_with_survival = scl, irreversible_mortality = irr,
                   shape = "varying", save_draws = TRUE)
      ce <- by_draw(do.call(predict_heat_injury, c(args, varying_method = "cumulative"))$draws)
      iv <- by_draw(do.call(predict_heat_injury, c(args, varying_method = "individual"))$draws)
      for (d in 1:3) {
        lab  <- paste("scaled", scl, "irreversible", irr, "draw", d)
        P    <- draw_P(s, d, tr$temp)
        ced  <- ce[ce$.draw == d, ]
        expect_gt(ced$hi[ced$time == 34], 0, label = lab)        # repair was partial
        expect_lt(ced$hi[ced$time == 34], ced$hi[ced$time == 12])
        expect_equal(ced$survival, ref_cumulative(tr$time, tr$temp, P, T_c = 26, rep_rate = rr,
                                                  scale = scl, irreversible = irr),
                     tolerance = 1e-9, info = lab)
        expect_equal(iv$survival[iv$.draw == d],
                     ref_individual(tr$time, tr$temp, P, T_c = 26, rep_rate = rr, scale = scl,
                                    irreversible = irr)$survival,
                     tolerance = 1e-10, info = lab)
      }
    }
  }
})

test_that("cumulative exposure: exact at a constant temperature and for a shared-injury mechanism", {
  # Every model, one temperature: survival is the fitted curve.
  spec_tail <- list(list(low = lin_T(rep(0.05, 3), 0), up = lin_T(rep(0.9, 3), 0),
                         k = exp_T(c(15.1, 15.4, 15.7), -0.68),
                         m0 = rep(1.4, 3), m1 = rep(-0.18, 3)))
  hot <- data.frame(time = 0:60, temp = 36)
  for (m in c("log-logistic", "log-normal", "weibull1", "weibull2")) {
    wf <- local_synthetic_wf(spec_tail, model = m)
    oc <- predict_heat_injury(hot, wf, shape = "varying", varying_method = "cumulative",
                              save_draws = TRUE)$draws
    P  <- draw_P(spec_tail[[1]], 1, hot$temp)
    curve_1 <- tdt_model_spec(m)$curve(P$low, P$up,
                                      P$k * log10(pmax(hot$time / 10^P$mid, 1e-12)))
    expect_equal(oc$survival[oc$.draw == 1], curve_1, tolerance = 1e-10, info = m)
  }
  # All individuals share one injury I, dI/dt = kappa r I^(1 - 1/kappa), and
  # thresholds are log-normal: at a constant T the curve is the log-normal model
  # (k = kappa ln 10, mid = -log10 r); under any trace the survival is exact.
  a0 <- 1.4; b0 <- -0.18; sa <- 0.2; sb <- 0.04
  kk <- function(temps) 1 / (sa + sb * (temps - 33))
  spec <- list(list(low = lin_T(c(0, 0), 0), up = lin_T(c(1, 1), 0),
                    k = function(temps) outer(c(1, 1), kk(temps)),
                    m0 = c(a0, a0), m1 = c(b0, b0)))
  wf <- local_synthetic_wf(spec, model = "log-normal")
  tm <- seq(0, 72, by = 0.25)
  tr <- data.frame(time = tm, temp = round(33 + 3 * sin(2 * pi * tm / 24), 4))
  kap <- kk(tr$temp) / log(10); r <- 10^(-(a0 + b0 * (tr$temp - 33)))
  I <- numeric(nrow(tr))
  for (j in 2:nrow(tr)) I[j] <- (I[j - 1]^(1 / kap[j - 1]) + r[j - 1] * diff(tm)[j - 1])^kap[j - 1]
  truth <- stats::pnorm(log(pmax(I, 1e-300)), lower.tail = FALSE)
  o <- predict_heat_injury(tr, wf, shape = "varying", varying_method = "cumulative",
                           save_draws = TRUE)$draws
  expect_equal(o$survival[o$.draw == 1], truth, tolerance = 1e-9)
})

test_that("both methods are the constant shape when the shape is flat (every model)", {
  flat <- list(list(low = lin_T(rep(0.05, 2), 0), up = lin_T(rep(0.9, 2), 0),
                    k = lin_T(c(6, 5.5), 0), m0 = c(1.4, 1.38), m1 = c(-0.18, -0.18)))
  tr <- data.frame(time = 0:96)
  tr$temp <- ifelse(tr$time %% 24 < 6, 36, 25)
  rp <- repair_pars_test(0.02)                 # repairs part of each day's injury
  for (m in c("log-logistic", "log-normal", "weibull1", "weibull2")) {
    wf <- local_synthetic_wf(flat, model = m)
    for (meth in c("individual", "cumulative")) {
      for (irr in c(TRUE, FALSE)) {
        for (tg in list("relative", 0.5)) {
          args <- list(tr, wf, target_surv = tg, T_c = 30, repair = TRUE, repair_pars = rp,
                       irreversible_mortality = irr, save_draws = TRUE)
          hc <- by_draw(do.call(predict_heat_injury, c(args, shape = "constant"))$draws)
          hv <- by_draw(do.call(predict_heat_injury,
                                c(args, shape = "varying", varying_method = meth))$draws)
          lab <- paste(m, meth, "irreversible", irr, "target", tg)
          expect_gt(max(hv$hi), 10, label = lab)               # injury accrues
          expect_lt(min(hv$survival), 0.85, label = lab)       # and kills
          expect_equal(hv$survival, hc$survival, tolerance = 1e-9, info = lab)
          expect_equal(hv$hi, hc$hi, tolerance = 1e-9, info = lab)
        }
      }
    }
  }
})

test_that("heat injury is the same under both methods unless repair scales with survival", {
  wf <- local_synthetic_wf(spec_moving)
  tm <- seq(0, 48, by = 0.25)
  tr <- data.frame(time = tm, temp = round(33 + 3 * sin(2 * pi * tm / 24), 4))
  run <- function(meth, ...) {
    by_draw(predict_heat_injury(tr, wf, target_surv = 0.6, shape = "varying",
                                varying_method = meth, save_draws = TRUE, ...)$draws)
  }
  i_d <- run("individual"); c_d <- run("cumulative")
  expect_equal(c_d$hi, i_d$hi, tolerance = 1e-12)
  expect_gt(max(abs(c_d$survival - i_d$survival)), 1e-3)
  rep_u <- list(repair = TRUE, repair_pars = repair_pars_test(0.02),
                repair_scales_with_survival = FALSE)
  i_u <- do.call(run, c("individual", rep_u)); c_u <- do.call(run, c("cumulative", rep_u))
  expect_equal(c_u$hi, i_u$hi, tolerance = 1e-12)
  rep_s <- utils::modifyList(rep_u, list(repair_scales_with_survival = TRUE))
  i_s <- do.call(run, c("individual", rep_s)); c_s <- do.call(run, c("cumulative", rep_s))
  expect_gt(max(abs(c_s$hi - i_s$hi)), 1e-3)     # each method's own survival scales it
})

test_that("the solver bisects when a Newton step does not halve the residual", {
  # Two bins whose curves share up, with the root on the second bin's low: Newton
  # alone stalls on the kink there (up to the iteration cap).
  P <- list(low = c(0.22908063675347973, 0.24184775652363896),
            up = c(0.24891559029346572, 0.24891559029346572),
            k = c(15.993181649483549, 19.636882558375074),
            mid = c(0.40723191175273832, -0.54793249815702449))
  res <- hi_rank_integrate(matrix(P$low), matrix(P$up), matrix(P$k), matrix(P$mid),
                           bin = c(1, 2, 2), dt = c(1, 1), below = rep(FALSE, 3),
                           unit_h = 1, mod = tdt_model_spec("log-normal"))
  expect_lt(res$max_iter, 100)
  want <- ref_individual(0:2, c(30, 31, 31),
                         data.frame(low = P$low[c(1, 2, 2)], up = P$up[c(1, 2, 2)],
                                    k = P$k[c(1, 2, 2)], mid = P$mid[c(1, 2, 2)]),
                         model = "log-normal")$survival
  expect_equal(res$survival[, 1], want, tolerance = 1e-12)
})

test_that("temperature order matters for cumulative exposure, not for the per-individual method", {
  wf <- local_synthetic_wf(spec_moving)
  hot_first  <- data.frame(time = (0:144) / 6, temp = c(rep(36, 72), rep(31, 73)))
  cool_first <- data.frame(time = (0:144) / 6, temp = c(rep(31, 72), rep(36, 73)))
  end_s <- function(tr, meth) {
    end_rows(predict_heat_injury(tr, wf, shape = "varying", varying_method = meth,
                                 save_draws = TRUE))$survival
  }
  expect_equal(end_s(hot_first, "individual"), end_s(cool_first, "individual"),
               tolerance = 1e-12)
  expect_gt(max(abs(end_s(hot_first, "cumulative") - end_s(cool_first, "cumulative"))), 1e-3)
})

test_that("at or below T_c survival does not change, even where the curve's top is lower", {
  spec <- list(list(low = lin_T(rep(0.05, 2), 0), up = lin_T(c(0.95, 0.96), -0.03),
                    k = lin_T(rep(6, 2), 0), m0 = c(1.4, 1.41), m1 = rep(-0.18, 2)))
  wf <- local_synthetic_wf(spec)                 # up(37) is lower than up(33)
  tr <- data.frame(time = 0:6, temp = c(33, 33, 37, 37, 33, 33, 33))
  for (meth in c("individual", "cumulative")) {
    d <- predict_heat_injury(tr, wf, T_c = 37, shape = "varying", varying_method = meth,
                             save_draws = TRUE)$draws
    d1 <- d[d$.draw == 1, ]
    expect_equal(d1$survival[4], d1$survival[3], tolerance = 1e-15, info = meth)
  }
})
