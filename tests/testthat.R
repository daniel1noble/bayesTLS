# Entry point for `testthat::test_dir()`.
#
# Run from the project root:
#   testthat::test_dir("tests/testthat")
#
# Set RUN_BRMS_TESTS=true in the environment to enable the brms-fitting
# integration tests; without that flag those tests are skipped so the suite
# runs in seconds.

library(testthat)
testthat::test_dir("tests/testthat")
