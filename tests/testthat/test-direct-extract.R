# Fast, brms-free tests touching the direct CTmax/z parameterisation that do NOT
# need a Stan fit.
#
# The former coefficient-reconstruction tests (silent-NA guard, C(T) oracle,
# absolute-fit threshold cells, per-draw z/CTmax) tested tdt_loglt /
# tdt_z_from_pars / tdt_ctmax_from_pars, which were RETIRED when all extraction
# moved onto the posterior_linpred engine (R/tls_engine.R). The C(T)/threshold
# arithmetic now lives in brms's `nlf(mid ~ ...)` + tls_eval_subpars(), so the
# C0/C_p conflation is structurally impossible rather than guarded. Coverage now:
#  - z/CTmax inversion math (fast): test-tdt-z-local-gate.R, test-tdt-ctmax-vectorised.R,
#    test-derive_z.R (all exercise the pure tls_local_z / tls_invert_logLT).
#  - direct-fit truth recovery, threshold modes, direct<->midpoint equivalence,
#    per-group/grouped output (gated): test-direct-fixture.R, test-unify-parity.R.

test_that("fit_4pl maps non-canonical duration_unit aliases consistently (finding 6)", {
  base <- data.frame(logd = log10(rep(c(1, 2, 4, 8), 3)),
                     temp_c = rep(c(-2, 0, 2), each = 4), n_surv = 5L, n_total = 10L)
  l10 <- function(u) {
    d <- base; attr(d, "tdt_meta") <- list(temp_mean = 35, duration_unit = u,
                                           response_type = "count")
    fit_4pl(d, ctmax = ~ 1, fit = FALSE)$meta$log10_tref
  }
  expect_equal(l10("h"),     l10("hours"))     # alias == canonical (was tm=1 fallback)
  expect_equal(l10("hr"),    l10("hours"))
  expect_equal(l10("Hours"), l10("hours"))     # case-insensitive
  expect_equal(l10("min"),   l10("minutes"))
})
