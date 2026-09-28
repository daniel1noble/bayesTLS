# Phylogenetic / covariance-matrix random effects (issue #12).
# Fast tier: formula parsing, gr()-aware grouping variables, data2 validation,
# spec construction, stancode smoke, and the re_formula auto-resolution rule.
# Plan + rationale: notes/2026-09-25-phylo-gr-support-plan.qmd.

ns <- function(f) getFromNamespace(f, "bayesTLS")

## Shared toy objects -----------------------------------------------------------
sp <- c("s1", "s2", "s3")
A  <- diag(3); A[A == 0] <- 0.4; dimnames(A) <- list(sp, sp)

make_std <- function() {
  raw <- expand.grid(species = sp, temp = c(38, 40, 42), t_hours = c(1, 4, 16))
  raw$n_total <- 20L
  set.seed(1)
  raw$n_surv  <- rbinom(nrow(raw), raw$n_total, 0.6)
  standardize_data(raw, temp = "temp", duration = "t_hours",
                   n_total = "n_total", n_surv = "n_surv",
                   duration_unit = "hours")
}

## T1: the AST term splitter ----------------------------------------------------
test_that("tdt_split_bars separates fixed terms from bar terms on every shape", {
  split <- ns("tdt_split_bars")

  s <- split("1 + (1 | gr(species, cov = A))")
  expect_identical(s$fixed, "1")
  expect_length(s$bars, 1L)

  s <- split("0 + region + (1 | gr(species, cov = A)) + (1 | batch)")
  expect_identical(s$fixed, "0 + region")
  expect_length(s$bars, 2L)

  s <- split("temp_c + (temp_c | site) + (1 | a/b)")
  expect_identical(s$fixed, "temp_c")
  expect_length(s$bars, 2L)

  s <- split("1 + (1 | Date) + (1 | Tank)")          # legacy plain terms
  expect_identical(s$fixed, "1")
  expect_length(s$bars, 2L)

  s <- split("temp_c * species")                      # no bars at all
  expect_identical(s$fixed, "temp_c * species")
  expect_length(s$bars, 0L)
})

test_that("rhs_fixed_vars survives gr() terms and keeps legacy behaviour", {
  rfv <- ns("rhs_fixed_vars")
  expect_identical(rfv("1 + (1 | gr(species, cov = A))"), character(0))
  expect_identical(rfv("0 + region + (1 | gr(species, cov = A))"), "region")
  expect_identical(rfv("temp_c + (temp_c | site)"), character(0))
  expect_identical(rfv("0 + species"), "species")          # legacy
  expect_identical(rfv("1 + (1 | Date) + (1 | G_Room)"), character(0))
})

## T2: grouping variables, gr()-aware ------------------------------------------
test_that("bar grouping variables exclude cov=/dist= and reject mm()", {
  bgv   <- ns("tdt_bar_group_vars")
  split <- ns("tdt_split_bars")
  bar1 <- split("(1 | gr(species, cov = A))")$bars[[1]]
  expect_identical(bgv(bar1), "species")                    # never "A"
  bar2 <- split("(1 | gr(g, by = b))")$bars[[1]]
  expect_identical(sort(bgv(bar2)), c("b", "g"))
  bar3 <- split("(1 | a/b)")$bars[[1]]
  expect_identical(sort(bgv(bar3)), c("a", "b"))
  bar4 <- split("(1 | mm(m1, m2))")$bars[[1]]
  expect_error(bgv(bar4), "mm\\(\\)")
})

## T3: normalize_cellmeans ------------------------------------------------------
test_that("normalize_cellmeans passes gr() formulas through unmangled and
           still rewrites single-factor treatment terms", {
  nc  <- ns("normalize_cellmeans")
  std <- make_std()

  f <- ~ 1 + (1 | gr(species, cov = A))
  expect_identical(nc(f, std), f)                           # untouched

  out <- nc(~ species + (1 | Date), std)                    # legacy rewrite
  rhs <- paste(deparse(out[[length(out)]]), collapse = " ")
  expect_match(rhs, "0 \\+ species")
  expect_match(rhs, "\\(1 \\| Date\\)")

  out2 <- nc(~ species + (1 | gr(species, cov = A)), std)   # rewrite + gr kept
  rhs2 <- paste(deparse(out2[[length(out2)]]), collapse = " ")
  expect_match(rhs2, "0 \\+ species")
  expect_match(rhs2, "gr\\(species, cov = A\\)")
})

## T5: data2 validation ---------------------------------------------------------
test_that("data2 validation errors are actionable", {
  val <- ns("tdt_validate_data2")
  std <- make_std()
  f   <- list(ctmax = ~ 1 + (1 | gr(species, cov = A)), z = ~ 1)

  expect_error(val(f, std, data2 = NULL), "data2")
  expect_error(val(f, std, data2 = list(B = A)), "A")
  bad <- A; bad[1, 2] <- 0.9
  expect_error(val(f, std, data2 = list(A = bad)), "symmetric")
  noname <- unname(A)
  expect_error(val(f, std, data2 = list(A = noname)), "dimnames|names")
  small <- A[1:2, 1:2]
  expect_error(val(f, std, data2 = list(A = small)), "level")
  cov_not_cor <- A * 2
  expect_warning(val(f, std, data2 = list(A = cov_not_cor)), "correlation")
  expect_silent(val(f, std, data2 = list(A = A)))
})

## T4 + T7: spec construction end to end (no Stan) ------------------------------
test_that("fit_4pl(fit = FALSE) builds a gr() spec; formula carries the term
           verbatim; missing pieces error clearly", {
  std <- make_std()
  f   <- ~ 1 + (1 | gr(species, cov = A))

  wf <- fit_4pl(std, ctmax = f, z = f, t_ref = 60,
                data2 = list(A = A), fit = FALSE)
  expect_s3_class(wf, "bayes_tls")
  ftxt <- paste(deparse(wf$formula), collapse = " ")
  expect_match(ftxt, "gr(species, cov = A)", fixed = TRUE)
  expect_identical(wf$data2, list(A = A))

  # gr() used but data2 not given -> actionable error
  expect_error(fit_4pl(std, ctmax = f, z = f, t_ref = 60, fit = FALSE),
               "data2")

  # grouping column absent from the data -> actionable error naming the column
  std2 <- std; std2$species <- NULL
  expect_error(fit_4pl(std2, ctmax = f, z = f, t_ref = 60,
                       data2 = list(A = A), fit = FALSE),
               "species")
})

## T6: brms accepts the built spec (codegen only, no compile) -------------------
test_that("make_stancode accepts the gr() spec and injects the covariance", {
  skip_if_not_installed("brms")
  std <- make_std()
  f   <- ~ 1 + (1 | gr(species, cov = A))
  wf  <- fit_4pl(std, ctmax = f, z = f, t_ref = 60,
                 data2 = list(A = A), fit = FALSE)
  # Muffle only brms's harmless redundant-catch-all note (the same one
  # fit_4pl() muffles at its brm() call); everything else must surface.
  code <- withCallingHandlers(
    brms::make_stancode(wf$formula, data = wf$data,
                        prior = wf$prior, data2 = list(A = A)),
    warning = function(w) {
      if (grepl("will not be used in the model", conditionMessage(w),
                fixed = TRUE))
        invokeRestart("muffleWarning")
    })
  expect_true(grepl("Lcov|cholesky", as.character(code), ignore.case = TRUE))
})

## T8: re_formula auto-resolution ----------------------------------------------
test_that("re_formula auto-rule: NA for fixed-effect by; the exact bar term for
           RE-only by; never blanket-include", {
  auto <- ns("tdt_auto_re_formula")

  pf <- list(ctmax = "1 + (1 | gr(species, cov = A)) + (1 | batch)",
             z     = "1 + (1 | gr(species, cov = A))")

  # by is a fixed-effect moderator elsewhere -> keep NA (legacy path)
  expect_identical(auto(pf, by = "species", fixed_vars = "species"), NA)

  # by exists only as a random effect -> exactly the matching bar term(s)
  rf <- auto(pf, by = "species", fixed_vars = character(0))
  expect_s3_class(rf, "formula")
  rtxt <- paste(deparse(rf), collapse = " ")
  expect_match(rtxt, "gr(species, cov = A)", fixed = TRUE)
  expect_false(grepl("batch", rtxt))                        # no blanket include

  # by matches nothing -> NA (engine will surface its own error downstream)
  expect_identical(auto(pf, by = "region", fixed_vars = character(0)), NA)
})
