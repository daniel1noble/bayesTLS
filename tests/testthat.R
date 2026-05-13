# Standard testthat entry point. R CMD check / devtools::test() /
# devtools::check() all run this file. test_check() auto-discovers
# test-*.R files inside tests/testthat/.
#
# To run from the project root in a plain R session:
#   testthat::test_dir("tests/testthat")
#
# Set RUN_BRMS_TESTS=true in the environment to enable the brms-fitting
# integration tests; without that flag those tests are skipped so the suite
# runs in seconds.

library(testthat)
library(bayesTLS)
test_check("bayesTLS")
