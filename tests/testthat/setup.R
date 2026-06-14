# Route any incidental plotting during the test run to a null device.
#
# Several tests exercise drawing code as a side effect — e.g. plot(<bayes_tls>)
# renders brms trace/density plots, and autoprinted ggplots draw to the active
# device. In a non-interactive R session with no device open, R falls back to
# writing "Rplots.pdf" into the working directory, which litters the working
# tree (tests/testthat/Rplots.pdf). Opening a null PDF device sends all of that
# output to nowhere. The device is closed automatically when the test session
# ends, so no explicit teardown (and no extra dependency) is needed.
grDevices::pdf(NULL)
