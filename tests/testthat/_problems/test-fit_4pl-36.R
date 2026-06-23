# Extracted from test-fit_4pl.R:36

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "bayesTLS", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
pf <- make_4pl_formula()$pforms
expect_match(paste(deparse(pf$lowraw), collapse = " "), "temp_c")
expect_match(paste(deparse(pf$upraw),  collapse = " "), "temp_c")
expect_match(paste(deparse(pf$logk),   collapse = " "), "temp_c")
expect_match(paste(deparse(pf$mid),    collapse = " "), "temp_c")
pf2 <- make_4pl_formula(temp_effects = "mid")$pforms
expect_equal(paste(deparse(pf2$lowraw), collapse = " "), "lowraw ~ 1")
expect_equal(paste(deparse(pf2$upraw),  collapse = " "), "upraw ~ 1")
expect_equal(paste(deparse(pf2$logk),   collapse = " "), "logk ~ 1")
expect_match(paste(deparse(pf2$mid),    collapse = " "), "mid ~ temp_c")
pf3 <- make_4pl_formula(temp_effects = "mid",
                          random_effects = "Date")$pforms
expect_match(paste(deparse(pf3$mid), collapse = " "), "\\(1 \\| Date\\)")
expect_error(make_4pl_formula(temp_effects = c("low", "up")),
               "mid.*must always carry")
