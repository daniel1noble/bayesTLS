# Fast formula-construction tests for the midpoint parameterisation's moderator
# handling (no brms). Mirror the situations in ms/case_studies_new.qmd:
#   (1) bare fit  -> species-agnostic, all four sub-parameters ~ temp_c
#   (2) explicit  -> low/up/k = ~ temp_c + species + species:temp_c (honoured)
#   (3) shorthand -> low/up/k = ~ temp_c * species (identical MODEL to (2))
#   (4) by=       -> ~ temp_c * species on ALL four sub-parameters
# Assert: the right RHS lands on each sub-parameter; the explicit expansion and
# the `*` shorthand build the SAME design; group_vars/grouped are recorded; and
# the "low/up/k vary but mid does not" footgun warns.

# A small standardised dataset carrying a `species` moderator.
std_sp <- local({
  raw <- expand.grid(temperature_C = c(34, 36, 38, 40),
                     exposure_h    = c(0.5, 1, 2, 4),
                     species       = c("A", "B", "C"))
  raw$n     <- 30L
  raw$alive <- 15L
  standardize_data(raw, temp = "temperature_C", duration = "exposure_h",
                   n_total = "n", n_surv = "alive", duration_unit = "hours")
})

# RHS of a sub-parameter's brms sub-formula, as a string.
sub_rhs <- function(f, par) paste(deparse(f$pforms[[par]][[3]]), collapse = " ")

# Numeric model-matrix columns for a sub-parameter, sorted by content so column
# NAME/ORDER differences (e.g. `temp_c:speciesB` vs `speciesB:temp_c`) don't
# count -- two formulas spanning the same design compare equal.
sub_design <- function(f, par, df) {
  rhs <- paste(deparse(f$pforms[[par]][[3]]), collapse = " ")   # drop the LHS
  m   <- stats::model.matrix(stats::as.formula(paste("~", rhs)), df)
  unname(m[, order(apply(m, 2, paste, collapse = "|")), drop = FALSE])
}

test_that("bare midpoint fit is species-agnostic (marginal 4PL, all four ~ temp_c)", {
  f <- make_4pl_formula()
  for (p in c("lowraw", "upraw", "logk", "mid"))
    expect_equal(sub_rhs(f, p), "temp_c")

  wf <- fit_4pl(std_sp, t_ref = 60, fit = FALSE)
  expect_equal(wf$meta$group_vars, character(0))
  expect_false(wf$meta$grouped)
})

test_that("explicit low/up/k moderator formulas are honoured (and recorded)", {
  expect_warning(
    f <- make_4pl_formula(low = ~ temp_c + species + species:temp_c,
                          up  = ~ temp_c + species + species:temp_c,
                          k   = ~ temp_c + species + species:temp_c),
    "mid does not.*POOLED")
  for (p in c("lowraw", "upraw", "logk"))
    expect_match(sub_rhs(f, p), "species")
  expect_equal(sub_rhs(f, "mid"), "temp_c")          # mid left pooled

  wf <- suppressWarnings(
    fit_4pl(std_sp, low = ~ temp_c + species + species:temp_c,
            up = ~ temp_c + species + species:temp_c,
            k = ~ temp_c + species + species:temp_c, t_ref = 60, fit = FALSE))
  expect_equal(wf$meta$group_vars, "species")
  expect_true(wf$meta$grouped)
})

test_that("`temp_c * species` shorthand builds the SAME model as the explicit expansion", {
  f_exp <- suppressWarnings(make_4pl_formula(
    low = ~ temp_c + species + species:temp_c,
    up  = ~ temp_c + species + species:temp_c,
    k   = ~ temp_c + species + species:temp_c))
  f_star <- suppressWarnings(make_4pl_formula(
    low = ~ temp_c * species, up = ~ temp_c * species, k = ~ temp_c * species))

  df <- data.frame(temp_c  = rep(c(-3, -1, 1, 3), 3),
                   species = factor(rep(c("A", "B", "C"), each = 4)))
  for (p in c("lowraw", "upraw", "logk"))
    expect_equal(sub_design(f_exp, p, df), sub_design(f_star, p, df))
})

test_that("by= applies temp_c * species to all four sub-parameters (overrides temp_effects)", {
  f <- make_4pl_formula(by = "species")
  for (p in c("lowraw", "upraw", "logk", "mid"))
    expect_equal(sub_rhs(f, p), "temp_c * species")

  # by also overrides a restricted temp_effects (it is the more specific request)
  f2 <- make_4pl_formula(temp_effects = "mid", by = "species")
  for (p in c("lowraw", "upraw", "logk", "mid"))
    expect_equal(sub_rhs(f2, p), "temp_c * species")

  wf <- fit_4pl(std_sp, by = "species", t_ref = 60, fit = FALSE)
  expect_equal(wf$meta$group_vars, "species")
  expect_true(wf$meta$grouped)
})

test_that("by= and the four explicit `~ temp_c * species` formulas build the same model", {
  f_by  <- make_4pl_formula(by = "species")
  f_exp <- make_4pl_formula(low = ~ temp_c * species, up = ~ temp_c * species,
                            k = ~ temp_c * species,  mid = ~ temp_c * species)
  df <- data.frame(temp_c  = rep(c(-3, -1, 1, 3), 3),
                   species = factor(rep(c("A", "B", "C"), each = 4)))
  for (p in c("lowraw", "upraw", "logk", "mid"))
    expect_equal(sub_design(f_by, p, df), sub_design(f_exp, p, df))
})
