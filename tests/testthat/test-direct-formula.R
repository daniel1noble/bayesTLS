# Tier-1 fast, deterministic, brms-free tests for the direct CTmax/z formula
# resolution. No model is fitted; everything is a formula-string assertion.

# RHS of a pform sub-formula (e.g. f$pforms$mid) as a string. Whitespace is
# collapsed because deparse() pads/wraps long expressions.
rhs_of <- function(f, par) {
  fm <- f$pforms[[par]]
  trimws(gsub("[[:space:]]+", " ", paste(deparse(fm[[length(fm)]]), collapse = " ")))
}
# All sub-formula strings (env-independent), for comparing two brmsformulas.
deparse_all <- function(f) {
  parts <- c(list(f$formula), f$pforms)
  sort(vapply(parts, function(x) paste(deparse(x), collapse = " "), character(1)))
}

test_that("midpoint mode is unchanged when neither ctmax nor z is supplied", {
  f <- make_4pl_formula()
  expect_setequal(names(f$pforms), c("lowraw", "upraw", "logk", "mid"))
  expect_false(any(c("CTmaxdev", "logz") %in% names(f$pforms)))
  # mid is the linear-in-T backbone, NOT reconstructed from CTmax/z
  expect_match(rhs_of(f, "mid"), "temp_c")
  expect_false(grepl("CTmaxdev", rhs_of(f, "mid")))
})

test_that("supplying ctmax/z triggers direct mode: CTmaxdev + logz, mid nlf-derived", {
  f <- make_4pl_formula(ctmax = ~ 1, z = ~ 1)
  expect_true(all(c("CTmaxdev", "logz") %in% names(f$pforms)))
  expect_match(rhs_of(f, "mid"), "CTmaxdev")        # mid reconstructed from CTmax & z
  expect_match(rhs_of(f, "mid"), "exp\\(logz\\)")
  expect_equal(rhs_of(f, "CTmaxdev"), "1")
  expect_equal(rhs_of(f, "logz"), "1")
  # intercept-only ctmax/z -> shape ~ temp_c
  expect_equal(rhs_of(f, "lowraw"), "temp_c")
  expect_equal(rhs_of(f, "logk"),   "temp_c")
})

test_that("cell-means FE is inherited onto the shape crossed with temp_c; RE is NOT", {
  f <- make_4pl_formula(ctmax = ~ 0 + life_stage + (1 | Date),
                        z     = ~ 0 + life_stage + (1 | Date))
  expect_equal(rhs_of(f, "lowraw"), "0 + life_stage + temp_c:life_stage")
  expect_equal(rhs_of(f, "upraw"),  "0 + life_stage + temp_c:life_stage")
  expect_equal(rhs_of(f, "logk"),   "0 + life_stage + temp_c:life_stage")
  # random effect stays on CTmax/z, never inherited onto shape
  for (sp in c("lowraw", "upraw", "logk")) expect_false(grepl("Date", rhs_of(f, sp)))
  expect_match(rhs_of(f, "CTmaxdev"), "\\(1 \\| Date\\)")
  expect_match(rhs_of(f, "logz"),     "\\(1 \\| Date\\)")
})

test_that("treatment-coded FE is inherited as G * temp_c", {
  f <- make_4pl_formula(ctmax = ~ sex)
  expect_equal(rhs_of(f, "lowraw"), "(sex) * temp_c")
})

test_that("an explicit shape formula overrides inheritance", {
  f <- make_4pl_formula(ctmax = ~ 0 + life_stage, low = ~ temp_c)
  expect_equal(rhs_of(f, "lowraw"), "temp_c")
  # the un-specified shapes still inherit
  expect_equal(rhs_of(f, "upraw"), "0 + life_stage + temp_c:life_stage")
})

test_that("pooled z: omitting z gives logz ~ 1 while ctmax keeps its structure", {
  f <- make_4pl_formula(ctmax = ~ 0 + sex)
  expect_equal(rhs_of(f, "logz"), "1")
  expect_equal(rhs_of(f, "CTmaxdev"), "0 + sex")
  # shape inherits from ctmax (the fallback source is z, but ctmax wins)
  expect_equal(rhs_of(f, "lowraw"), "0 + sex + temp_c:sex")
})

test_that("threshold='absolute' adds the C(T) correction to mid; relative omits it", {
  fr <- make_4pl_formula(ctmax = ~ 1, threshold = "relative")
  fa <- make_4pl_formula(ctmax = ~ 1, threshold = "absolute")
  expect_false(grepl("log\\(\\(up", rhs_of(fr, "mid")))
  expect_match(rhs_of(fa, "mid"), "log\\(\\(up - 0.5\\)/\\(0.5 - low\\)\\)")
  # non-default p flows into the correction
  fp <- make_4pl_formula(ctmax = ~ 1, threshold = "absolute", p = 0.9)
  expect_match(rhs_of(fp, "mid"), "log\\(\\(up - 0.9\\)/\\(0.9 - low\\)\\)")
})

test_that("resolve_shape() truth table (named internal helper, D-T3)", {
  expect_equal(resolve_shape(NULL, "1"), "temp_c")
  expect_equal(resolve_shape(NULL, "0 + life_stage + (1 | Date)"),
               "0 + life_stage + temp_c:life_stage")
  expect_equal(resolve_shape(NULL, "0 + life_stage"),
               "0 + life_stage + temp_c:life_stage")
  expect_equal(resolve_shape(NULL, "sex"), "(sex) * temp_c")
  expect_equal(resolve_shape(NULL, "1"), "temp_c")
  expect_equal(resolve_shape(~ temp_c, "0 + life_stage"), "temp_c")  # explicit wins
})

test_that("bounds = c(lower, upper) supersedes lower/upper with identical constants", {
  expect_equal(deparse_all(make_4pl_formula(bounds = c(0.85, 1))),
               deparse_all(make_4pl_formula(lower = 0.85, upper = 1)))
  # the sublethal lower bound is baked into the asymptote constant
  f <- make_4pl_formula(bounds = c(0.85, 1))
  expect_match(paste(deparse(f$formula), collapse = " "), "0.85")
})

test_that("a non-formula passed to a formula arg errors clearly", {
  expect_error(make_4pl_formula(ctmax = 0.5), "one-sided formula")
})
