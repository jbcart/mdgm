test_that("standalone model creation works", {
  edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree")
  expect_s3_class(model, "SrfModel")
  expect_equal(model$nvertices(), 3L)
  expect_equal(model$ncolors(), 2L)
  expect_false(model$has_emission())
})

test_that("hierarchical model creation works", {
  edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree",
                       emission = "bernoulli")
  expect_true(model$has_emission())
})

test_that("acyclic orientation model works", {
  edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "acyclic_orientation")
  expect_equal(model$nvertices(), 3L)
})

test_that("invalid dag_type is rejected", {
  edges <- rbind(c(1, 2), c(2, 1))
  nug <- nug_from_edge_list(2, edges, seed = 42L)
  expect_error(mdgm_model(nug, dag_type = "invalid"))
})

test_that("MdgmModel alias still works", {
  edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree")
  expect_s3_class(model, "SrfModel")
  # MdgmModel is an alias for SrfModel
  expect_true(inherits(model, "SrfModel"))
})

test_that("MRF model creation works (exchange)", {
  edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- srf_model(nug, spatial = mrf(method = "exchange"))
  expect_s3_class(model, "SrfModel")
  expect_equal(model$model_type(), "mrf")
  expect_equal(model$nvertices(), 3L)
  expect_equal(model$ncolors(), 2L)
  expect_false(model$has_emission())
})

test_that("MRF model creation works (pseudo_likelihood + emission)", {
  edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- srf_model(nug, spatial = mrf(method = "pseudo_likelihood"),
                     emission = "bernoulli", n_colors = 3L)
  expect_true(model$has_emission())
  expect_equal(model$ncolors(), 3L)
  expect_equal(model$emission_type(), "bernoulli")
})

test_that("invalid spatial config is rejected", {
  edges <- rbind(c(1, 2), c(2, 1))
  nug <- nug_from_edge_list(2, edges, seed = 42L)
  expect_error(srf_model(nug, spatial = list(foo = "bar")),
               "spatial must be created by mdgm")
})
