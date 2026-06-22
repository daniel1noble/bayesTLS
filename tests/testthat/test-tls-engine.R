# Fast, brms-free tests for the engine plumbing in R/tls_engine.R that does not
# call posterior_linpred: the grid builder and the by-resolution. (tls_local_z /
# tls_invert_logLT are covered in test-tdt-z-local-gate.R and
# test-tdt-ctmax-vectorised.R; tls_eval_subpars / tls_draw_ids are exercised
# gated through every migrated reader on the cached fixtures.)

test_that("tls_build_grid: ungrouped grid is temp_grid long, fills model columns, .grp = 'all'", {
  mdata <- data.frame(temp_c = c(-2, 0, 2), logd = 1, life_stage = factor("a"))
  nd <- tls_build_grid(mdata, by = NULL, temp = "temp_c", temp_grid = c(-1, 0, 1))
  expect_equal(nd$temp_c, c(-1, 0, 1))
  expect_true(all(nd$.grp == "all"))
  expect_true("life_stage" %in% names(nd))           # missing model column filled from row 1
  expect_true(all(nd$life_stage == "a"))
})

test_that("tls_build_grid: grouped grid crosses moderator levels x temp_grid, tagged in .grp", {
  mdata <- data.frame(temp_c = c(-2, 0, 2),
                      life_stage = factor(c("young", "old", "young")))
  nd <- tls_build_grid(mdata, by = "life_stage", temp = "temp_c", temp_grid = c(0, 1))
  expect_setequal(unique(nd$.grp), c("young", "old"))
  expect_equal(nrow(nd), 2L * 2L)                     # 2 levels x 2 temps
  expect_equal(sum(nd$life_stage == "young"), 2L)
  expect_setequal(unique(nd$temp_c), c(0, 1))
})

test_that("tdt_resolve_by: explicit wins; NULL uses meta$group_vars; single-condition -> NULL", {
  wf_g <- structure(list(meta = list(group_vars = "sex")), class = "bayes_tls")
  wf_s <- structure(list(meta = list(group_vars = character(0))), class = "bayes_tls")
  expect_equal(tdt_resolve_by(wf_g, by = NULL), "sex")       # auto from the fit
  expect_equal(tdt_resolve_by(wf_g, by = "clone"), "clone")  # explicit overrides
  expect_null(tdt_resolve_by(wf_s, by = NULL))               # single-condition
})
