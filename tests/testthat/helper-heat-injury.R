# Reference for predict_heat_injury(shape = "varying"), written out from the
# help page for ONE parameter set, independently of the package integrator.
#
# An individual of survival rank p dies, at a constant temperature T, when
# population survival there falls to p: at LT_p(T) = 10^(mid(T) + c_p(T)), c_p
# the model's absolute offset. Injury is additive per individual,
# sum(dt / LT_p(T)), and the individual dies when it reaches 1. Outstanding
# relative dose is kept by the temperature at which it was accrued (A_T), so
# individual p's injury is sum_T A_T LT_rel(T) / LT_p(T); repair removes the
# same fraction from every A_T. Survival is the rank whose injury is exactly 1
# (running minimum when mortality is irreversible); heat injury is the relative
# dose, or for an absolute target the injury of that rank. Each interval is
# held at its starting temperature. Survival starts at the top of the curve at
# the first point whose interval accrues injury; ranks at or above a visited
# temperature's up die there for good (reversible mortality rebounds only up to
# the lowest such up); repair scaled by survival divides by the starting up.
# `P` holds low, up, k, mid at each trace point; found by uniroot, so it shares
# no code with the package's solver.
ref_individual <- function(time, temp, P, model = "log-logistic", p_target = NULL,
                           conv = 1, T_c = NULL, rep_rate = NULL, scale = TRUE,
                           irreversible = TRUE) {
  spec <- tdt_model_spec(model)
  n <- length(time); dt <- diff(time)
  ut <- unique(temp); bin <- match(temp, ut)
  Pb <- P[match(ut, temp), , drop = FALSE]        # the curve at each bin's temperature
  below <- if (is.null(T_c)) rep(FALSE, n) else temp <= T_c
  A <- numeric(length(ut))
  injury <- function(p) {
    act <- A > 0
    if (any(act & p >= Pb$up)) return(Inf)        # passed: threshold time 0
    ok <- act & Pb$low < p
    sum(A[ok] * 10^(-spec$abs_corr(Pb$low[ok], Pb$up[ok], Pb$k[ok], p)))
  }
  rr <- if (is.null(rep_rate)) numeric(n) else rep_rate
  i1 <- which(dt > 0 & !below[-n]); i1 <- if (length(i1)) i1[1] else 1
  S <- H <- numeric(n)
  S[1] <- spec$curve(P$low[i1], P$up[i1], P$k[i1] * log10(1e-12))
  cap <- S[1]
  for (j in 2:n) {
    a <- j - 1; w <- dt[a]; b <- bin[a]
    if (w > 0 && !below[a]) {
      A[b] <- A[b] + w / 10^Pb$mid[b]
      cap <- min(cap, Pb$up[b])
    }
    if (!is.null(rep_rate) && w > 0) {
      R <- rr[a] * w * conv * (if (scale) S[a] / P$up[i1] else 1)
      tot <- sum(A)
      if (tot > 0) A <- if (R >= tot) 0 * A else A * (1 - R / tot)
    }
    act <- A > 0
    root <- 1
    if (any(act)) {
      lo <- min(Pb$low[act]); hi <- min(Pb$up[act])
      root <- hi
      if (lo < hi) {
        top <- hi - 1e-13 * (hi - lo)
        f <- function(p) injury(p) - 1
        root <- if (f(top) < 0) top else stats::uniroot(f, c(lo, top), tol = 1e-15)$root
      }
    }
    S[j] <- if (irreversible) min(S[a], root) else min(cap, root)
    H[j] <- if (is.null(p_target)) sum(A) else injury(p_target)
  }
  data.frame(dose = H, survival = S)
}

# Reference for varying_method = "cumulative" (cumulative exposure), for ONE
# parameter set, written from its definition: survival is the state. A step
# that injures moves the state to that step's curve at the relative dose where
# the curve has the injury's survival (zero dose at or above the curve's top),
# unless that survival is at or below the curve's floor: then the state stays
# on its own curve and the step adds no injury to it. Repair lowers the dose on
# the curve the state lies on; at or below T_c nothing moves and only repair
# acts. Survival is the injury's survival, as a running minimum (irreversible)
# or capped at the lowest up of the temperatures that injured (reversible).
# Survival starts at the top of the curve at the first point whose interval
# accrues injury, and repair scaled by survival divides by that up. Doses come
# from survival through the model's offset (abs_corr), so this shares no code
# with the package's mapping; use traces where survival moves measurably away
# from the curve's top.
ref_cumulative <- function(time, temp, P, model = "log-logistic", conv = 1,
                           T_c = NULL, rep_rate = NULL, scale = TRUE,
                           irreversible = TRUE) {
  spec <- tdt_model_spec(model)
  n <- length(time); dt <- diff(time)
  below <- if (is.null(T_c)) rep(FALSE, n) else temp <= T_c
  surv_at <- function(D, j) spec$curve(P$low[j], P$up[j], P$k[j] * log10(max(D, 1e-12)))
  dose_at <- function(s, j) {                     # the dose at which curve j has survival s
    if (s >= P$up[j]) return(0)
    10^spec$abs_corr(P$low[j], P$up[j], P$k[j], s)
  }
  rr <- if (is.null(rep_rate)) numeric(n) else rep_rate
  i1 <- which(dt > 0 & !below[-n]); i1 <- if (length(i1)) i1[1] else 1
  S <- numeric(n); S[1] <- surv_at(0, i1); cap <- S[1]
  cb <- i1; Dc <- 0                               # the state's curve (a trace point) and dose
  for (j in 2:n) {
    a <- j - 1; w <- dt[a]
    if (w <= 0) { S[j] <- S[a]; next }
    on <- FALSE
    if (!below[a]) {
      cap <- min(cap, P$up[a])
      if (temp[a] == temp[cb]) {
        on <- TRUE
      } else {
        sc <- surv_at(Dc, cb)                     # the injury's survival
        if (sc > P$low[a]) { Dc <- dose_at(sc, a); cb <- a; on <- TRUE }
      }
    }
    R  <- rr[a] * w * conv * (if (scale) S[a] / P$up[i1] else 1)
    Dc <- max(0, Dc + (if (on) w / 10^P$mid[a] else 0) - R)
    S[j] <- min(if (irreversible) S[a] else cap, surv_at(Dc, cb))
  }
  S
}
