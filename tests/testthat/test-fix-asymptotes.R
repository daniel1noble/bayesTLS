# fit_4pl(fix_asymptotes = ) / make_4pl_priors(fix_asymptotes = ): pin one or
# both asymptotes through a brms constant() prior with the sub-parameter
# formula forced to `~ 1`. Default (NULL) must be byte-identical to before.
# Fast tests need no Stan; the integration tests at the end are gated.

fa_std <- function(seed = 1) {
  raw <- simulate_tdt(seed = seed)
  standardize_data(raw, temp = "T", duration = "t_hours",
                   n_total = "n_trials", n_surv = "y_alive",
                   duration_unit = "hours")
}
fa_std2 <- function(seed = 2) {
  raw <- simulate_tdt_2group(seed = seed)
  standardize_data(raw, temp = "T", duration = "t_hours",
                   n_total = "n_trials", n_surv = "y_alive",
                   duration_unit = "hours")
}
pf_rhs <- function(wf, nlpar) {
  paste(deparse(wf$formula$pforms[[nlpar]][[3]]), collapse = "")
}
asy_rows <- function(pr, nlpar) {
  d <- as.data.frame(pr)
  d <- d[d$nlpar == nlpar, c("prior", "class", "coef", "nlpar")]
  rownames(d) <- NULL
  d
}
constant_value <- function(prior_string) {
  as.numeric(sub("^constant\\((.*)\\)$", "\\1", prior_string))
}
is_constant_intercept <- function(rows) {
  nrow(rows) == 1L && rows$coef == "Intercept" &&
    grepl("^constant\\(", rows$prior)
}

# ---- resolver and raw mapping ------------------------------------------------

test_that("tdt_resolve_fix_asymptotes accepts the documented forms", {
  expect_null(tdt_resolve_fix_asymptotes(NULL))
  expect_null(tdt_resolve_fix_asymptotes(FALSE))
  expect_equal(tdt_resolve_fix_asymptotes(TRUE), list(low = 0, up = 1))
  expect_equal(tdt_resolve_fix_asymptotes(TRUE, 0.85, 1),
               list(low = 0.85, up = 1))
  expect_equal(tdt_resolve_fix_asymptotes("up"), list(up = 1))
  expect_equal(tdt_resolve_fix_asymptotes(c("low", "up")),
               list(low = 0, up = 1))
  expect_equal(tdt_resolve_fix_asymptotes(c(low = 0.02)), list(low = 0.02))
  expect_equal(tdt_resolve_fix_asymptotes(c(up = 0.9, low = 0)),
               list(up = 0.9, low = 0))
})

test_that("tdt_resolve_fix_asymptotes rejects malformed input", {
  expect_error(tdt_resolve_fix_asymptotes("mid"), "fix_asymptotes")
  expect_error(tdt_resolve_fix_asymptotes(c(mid = 1)), "fix_asymptotes")
  expect_error(tdt_resolve_fix_asymptotes(c(0, 1)), "fix_asymptotes") # unnamed
  expect_error(tdt_resolve_fix_asymptotes(c(low = 0, low = 0)),
               "fix_asymptotes")
  expect_error(tdt_resolve_fix_asymptotes(c(up = 1.2)),
               "outside the response bounds")
  expect_error(tdt_resolve_fix_asymptotes(c(low = -0.1)),
               "outside the response bounds")
  expect_error(tdt_resolve_fix_asymptotes(c(low = NA_real_)),
               "outside the response bounds")
  expect_error(tdt_resolve_fix_asymptotes(list(low = 0)), "fix_asymptotes")
})

test_that("tdt_fixed_asymptote_raw maps interior values exactly, pins bounds", {
  b <- compute_4pl_bounds(0, 1)
  for (v in c(0.02, 0.1, 0.3)) {
    r <- tdt_fixed_asymptote_raw(v, "low", b)
    expect_equal(r$value, v, tolerance = 1e-10)
    expect_equal(b$low_min + stats::plogis(r$raw) * b$low_w, v,
                 tolerance = 1e-10)
  }
  for (v in c(0.7, 0.9, 0.98)) {
    r <- tdt_fixed_asymptote_raw(v, "up", b)
    expect_equal(r$value, v, tolerance = 1e-10)
  }
  # the response bounds themselves sit at the padded window edge
  expect_lt(abs(tdt_fixed_asymptote_raw(0, "low", b)$value - b$low_min), 1e-8)
  expect_lt(abs(tdt_fixed_asymptote_raw(1, "up",  b)$value - b$up_max),  1e-8)
  expect_lt(abs(tdt_fixed_asymptote_raw(0, "low", b)$value - 0.001), 1e-8)
  expect_lt(abs(tdt_fixed_asymptote_raw(1, "up",  b)$value - 0.999), 1e-8)
  # sublethal bounds: the bound maps to that window's edge
  b2 <- compute_4pl_bounds(0.85, 1)
  r85 <- tdt_fixed_asymptote_raw(0.85, "low", b2, 0.85, 1)
  expect_lt(abs(r85$value - b2$low_min), 1e-8)
  # a value across the disjoint-bounds midpoint is an error, not a silent clamp
  expect_error(tdt_fixed_asymptote_raw(0.7, "low", b), "outside the interval")
  expect_error(tdt_fixed_asymptote_raw(0.3, "up",  b), "outside the interval")
})

# ---- default unchanged -------------------------------------------------------

test_that("the default (no fixing) is unchanged", {
  std <- fa_std()
  p0 <- make_4pl_priors(std)
  expect_identical(as.data.frame(make_4pl_priors(std, fix_asymptotes = NULL)),
                   as.data.frame(p0))
  expect_identical(as.data.frame(make_4pl_priors(std, fix_asymptotes = FALSE)),
                   as.data.frame(p0))
  wf <- fit_4pl(std, fit = FALSE)
  expect_null(wf$meta$fix_asymptotes)
  expect_identical(as.data.frame(wf$prior), as.data.frame(p0))
  expect_equal(pf_rhs(wf, "lowraw"), "temp_c")
  expect_equal(pf_rhs(wf, "upraw"),  "temp_c")
  expect_true(all(c("low", "up") %in% wf$meta$temp_effects))
})

# ---- fixing both -------------------------------------------------------------

test_that("fix_asymptotes = TRUE pins low and up at the padded bounds", {
  std <- fa_std()
  wf  <- fit_4pl(std, fix_asymptotes = TRUE, fit = FALSE)
  expect_equal(pf_rhs(wf, "lowraw"), "1")
  expect_equal(pf_rhs(wf, "upraw"),  "1")
  expect_equal(pf_rhs(wf, "logk"),   "temp_c")
  expect_equal(pf_rhs(wf, "mid"),    "temp_c")
  expect_true(is_constant_intercept(asy_rows(wf$prior, "lowraw")))
  expect_true(is_constant_intercept(asy_rows(wf$prior, "upraw")))
  fixed <- unname(wf$meta$fix_asymptotes[c("low", "up")])
  expect_lt(max(abs(fixed - c(0.001, 0.999))), 1e-8)
  expect_setequal(wf$meta$temp_effects, c("k", "mid"))
  # the constant reproduces the pinned value through the reparameterisation
  b   <- wf$meta$bounds
  raw <- constant_value(asy_rows(wf$prior, "upraw")$prior)
  expect_lt(abs(b$up_min + stats::plogis(raw) * b$up_w - 0.999), 1e-8)
  # everything else in the prior table is the default
  p0  <- as.data.frame(make_4pl_priors(std))
  keep <- function(d) {
    d <- d[!(d$nlpar %in% c("lowraw", "upraw")), ]
    rownames(d) <- NULL
    d
  }
  expect_identical(keep(as.data.frame(wf$prior)), keep(p0))
  expect_output(print(wf), "Fixed:.*low = 0.001, up = 0.999")
})

# ---- fixing ONE: the other must stay estimated ------------------------------

test_that("fixing up leaves low estimated (priors, formula, meta)", {
  std <- fa_std()
  p0  <- make_4pl_priors(std)
  for (spec in list(c(up = 1), "up")) {
    wf <- fit_4pl(std, fix_asymptotes = spec, fit = FALSE)
    expect_true(is_constant_intercept(asy_rows(wf$prior, "upraw")))
    expect_identical(asy_rows(wf$prior, "lowraw"), asy_rows(p0, "lowraw"))
    expect_equal(pf_rhs(wf, "lowraw"), "temp_c")   # still varies with T
    expect_equal(pf_rhs(wf, "upraw"),  "1")
    expect_equal(names(wf$meta$fix_asymptotes), "up")
    expect_true("low" %in% wf$meta$temp_effects)
    expect_false("up" %in% wf$meta$temp_effects)
  }
})

test_that("fixing low leaves up estimated (priors, formula, meta)", {
  std <- fa_std()
  p0  <- make_4pl_priors(std)
  for (spec in list(c(low = 0), "low")) {
    wf <- fit_4pl(std, fix_asymptotes = spec, fit = FALSE)
    expect_true(is_constant_intercept(asy_rows(wf$prior, "lowraw")))
    expect_identical(asy_rows(wf$prior, "upraw"), asy_rows(p0, "upraw"))
    expect_equal(pf_rhs(wf, "upraw"),  "temp_c")
    expect_equal(pf_rhs(wf, "lowraw"), "1")
    expect_equal(names(wf$meta$fix_asymptotes), "low")
    expect_true("up" %in% wf$meta$temp_effects)
    expect_false("low" %in% wf$meta$temp_effects)
  }
})

test_that("named values pin at the requested value", {
  std <- fa_std()
  wf  <- fit_4pl(std, fix_asymptotes = c(low = 0.02, up = 0.9), fit = FALSE)
  expect_equal(wf$meta$fix_asymptotes, c(low = 0.02, up = 0.9),
               tolerance = 1e-9)
  b   <- wf$meta$bounds
  raw <- constant_value(asy_rows(wf$prior, "lowraw")$prior)
  expect_equal(b$low_min + stats::plogis(raw) * b$low_w, 0.02, tolerance = 1e-9)
})

# ---- direct parameterisation and `by` ----------------------------------------

test_that("direct cell-means: fixed one is a constant, other stays per group", {
  std2 <- fa_std2()
  wf <- fit_4pl(std2, ctmax = ~ 0 + grp, z = ~ 0 + grp,
                fix_asymptotes = c(up = 1), fit = FALSE)
  expect_true(is_constant_intercept(asy_rows(wf$prior, "upraw")))
  low_rows <- asy_rows(wf$prior, "lowraw")
  expect_true(all(c("grpA", "grpB") %in% low_rows$coef)) # low still per group
  expect_equal(pf_rhs(wf, "upraw"), "1")
  expect_match(pf_rhs(wf, "lowraw"), "grp")
  expect_match(pf_rhs(wf, "logk"), "grp")
  wf2 <- fit_4pl(std2, ctmax = ~ 0 + grp, z = ~ 0 + grp,
                 fix_asymptotes = c(low = 0), fit = FALSE)
  expect_true(is_constant_intercept(asy_rows(wf2$prior, "lowraw")))
  expect_true(all(c("grpA", "grpB") %in% asy_rows(wf2$prior, "upraw")$coef))
})

test_that("by = moderator: only the fixed asymptote loses the by structure", {
  std2 <- fa_std2()
  wf <- fit_4pl(std2, by = "grp", fix_asymptotes = c(low = 0), fit = FALSE)
  expect_equal(pf_rhs(wf, "lowraw"), "1")
  expect_match(pf_rhs(wf, "upraw"), "temp_c \\* grp")
  expect_match(pf_rhs(wf, "mid"),   "temp_c \\* grp")
})

test_that("a conflicting explicit formula is refused", {
  std <- fa_std()
  expect_error(
    fit_4pl(std, fix_asymptotes = c(low = 0), low = ~ temp_c, fit = FALSE),
    "fixed by `fix_asymptotes` but also has the formula")
  expect_error(
    fit_4pl(std, fix_asymptotes = TRUE, up = ~ 1 + temp_c, fit = FALSE),
    "fixed by `fix_asymptotes`")
  # an explicit `~ 1` agrees with fixing and is accepted
  expect_s3_class(fit_4pl(std, fix_asymptotes = TRUE, low = ~ 1, fit = FALSE),
                  "bayes_tls")
  expect_error(fit_4pl(std, fix_asymptotes = c(low = 0.7), fit = FALSE),
               "outside the interval")
})

test_that("Stan code pins only the fixed coefficient", {
  skip_if_not_installed("brms")
  std <- fa_std()
  code_for <- function(spec) {
    wf <- fit_4pl(std, fix_asymptotes = spec, fit = FALSE)
    brms::make_stancode(wf$formula, data = std, prior = wf$prior,
                        family = brms::beta_binomial(link = "identity"))
  }
  both <- code_for(TRUE)
  expect_match(both, "b_lowraw\\[1\\] = -20\\.72")
  expect_match(both, "b_upraw\\[1\\] = 20\\.72")
  up_only <- code_for(c(up = 1))
  expect_match(up_only, "b_upraw\\[1\\] = 20\\.72")
  expect_false(grepl("b_lowraw\\[1\\] =", up_only))
  low_only <- code_for("low")
  expect_match(low_only, "b_lowraw\\[1\\] = -20\\.72")
  expect_false(grepl("b_upraw\\[1\\] =", low_only))
  # the mean function itself is unchanged by fixing
  default <- code_for(NULL)
  mu_line <- function(code)
    grep("mu\\[n\\] = ", strsplit(code, "\n")[[1]], value = TRUE)
  expect_identical(mu_line(both), mu_line(default))
})

# ---- gated integration: real fits --------------------------------------------

fa_fit <- function(std, ...) {
  fit_4pl(std, ..., chains = 2, iter = 800, warmup = 400, cores = 2, seed = 1,
          control = list(adapt_delta = 0.95, max_treedepth = 12), refresh = 0)
}
raw_draws <- function(wf, par) {
  as.numeric(posterior::as_draws_matrix(get_brmsfit(wf), variable = par))
}
# natural-scale asymptote draws from the raw Intercept draws
asy_draws <- function(wf, which) {
  b <- wf$meta$bounds
  if (which == "low")
    b$low_min + stats::plogis(raw_draws(wf, "b_lowraw_Intercept")) * b$low_w
  else
    b$up_min + stats::plogis(raw_draws(wf, "b_upraw_Intercept")) * b$up_w
}

test_that("fitted, up fixed: up constant over draws, low estimated, tls ok", {
  skip_unless_brms()
  std <- fa_std()
  wf  <- fa_fit(std, fix_asymptotes = c(up = 1))
  up_draws  <- asy_draws(wf, "up")
  low_draws <- asy_draws(wf, "low")
  expect_equal(stats::sd(up_draws), 0, tolerance = 1e-12)
  expect_equal(unique(round(up_draws, 6)), 0.999)
  expect_gt(stats::sd(low_draws), 0)
  expect_gt(stats::sd(raw_draws(wf, "b_lowraw_temp_c")), 0) # T slope estimated
  res <- tls(wf, params = c("z", "ctmax"))$summary
  expect_gte(nrow(res), 2L)
  expect_true(all(is.finite(res$median)))
  abs_res <- tls(wf, params = c("z", "ctmax"), target_surv = "absolute")$summary
  expect_true(all(is.finite(abs_res$median)))
})

test_that("fitted, low fixed: low constant across draws, up estimated", {
  skip_unless_brms()
  std <- fa_std()
  wf  <- fa_fit(std, fix_asymptotes = "low")
  low_draws <- asy_draws(wf, "low")
  up_draws  <- asy_draws(wf, "up")
  expect_equal(stats::sd(low_draws), 0, tolerance = 1e-12)
  expect_equal(unique(round(low_draws, 6)), 0.001)
  expect_gt(stats::sd(up_draws), 0)
  expect_output(print(wf), "Fixed:.*low = 0.001")
  res <- tls(wf, params = c("z", "ctmax"))$summary
  expect_true(all(is.finite(res$median)))
})
