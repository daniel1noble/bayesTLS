# Route any incidental plotting during the test run to a null device.
#
# Several tests exercise drawing code as a side effect — e.g. plot(<bayes_tls>)
# renders brms trace/density plots, and autoprinted ggplots draw to the active
# device. In a non-interactive R session with no device open, R falls back to
# writing "Rplots.pdf" into the working directory, which litters the working
# tree (tests/testthat/Rplots.pdf). Opening a null PDF device sends all of that
# output to nowhere. Closed again at the end of the test session.
grDevices::pdf(NULL)
if (requireNamespace("withr", quietly = TRUE)) {
  withr::defer(grDevices::dev.off(), teardown_env())
}
