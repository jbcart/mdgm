test_that("standalone MCMC runs and returns correct structure", {
  edges <- rbind(
    c(1, 2), c(2, 1), c(2, 3), c(3, 2),
    c(1, 3), c(3, 1), c(3, 4), c(4, 3)
  )
  nug <- nug_from_edge_list(4, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree")

  result <- mcmc(model,
                 z_init = c(0L, 0L, 1L, 1L),
                 psi_init = 0.5,
                 n_iter = 50L,
                 psi_tune = 0.1,
                 seed = 123L)

  expect_s3_class(result, "MdgmResult")

  # Check z dimensions
  z_mat <- result$z()
  expect_equal(dim(z_mat), c(4L, 50L))

  # Check psi length
  expect_length(result$psi(), 50L)
  expect_true(all(is.finite(result$psi())))

  # Check single iteration access
  z1 <- result$z(1)
  expect_length(z1, 4L)

  # Check DAG storage
  dag_mat <- result$dag()
  expect_equal(dim(dag_mat), c(4L, 50L))
})

test_that("hierarchical MCMC runs with Bernoulli emission", {
  edges <- rbind(
    c(1, 2), c(2, 1), c(2, 3), c(3, 2),
    c(1, 3), c(3, 1), c(3, 4), c(4, 3)
  )
  nug <- nug_from_edge_list(4, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree",
                       emission = "bernoulli")

  y <- list(c(1L, 0L, 1L), c(0L, 0L), c(1L, 1L, 1L), c(0L, 1L))

  result <- mcmc(model, y = y,
                 z_init = c(0L, 0L, 1L, 1L),
                 psi_init = 0.5,
                 eta_init = c(0.3, 0.7),
                 n_iter = 50L,
                 psi_tune = 0.1,
                 emission_prior_params = c(1, 1),
                 seed = 123L)

  expect_s3_class(result, "MdgmResult")

  # Check eta dimensions
  eta_mat <- result$eta()
  expect_equal(dim(eta_mat), c(2L, 50L))

  # Eta values should be between 0 and 1
  expect_true(all(eta_mat >= 0 & eta_mat <= 1))
})

test_that("acceptance rates are computed", {
  edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree")

  result <- mcmc(model, z_init = c(0L, 0L, 1L),
                 psi_init = 0.5, n_iter = 100L,
                 psi_tune = 0.1, seed = 42L)

  rates <- result$acceptance_rates()
  expect_named(rates, c("psi", "graph"))
  expect_true(rates["psi"] >= 0 && rates["psi"] <= 1)
})

test_that("summary prints without error", {
  edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree")

  result <- mcmc(model, z_init = c(0L, 0L, 1L),
                 psi_init = 0.5, n_iter = 50L,
                 psi_tune = 0.1, seed = 42L)

  expect_output(result$summary(), "MDGM MCMC Results")
})

test_that("edge inclusion probabilities work", {
  edges <- rbind(
    c(1, 2), c(2, 1), c(2, 3), c(3, 2),
    c(1, 3), c(3, 1)
  )
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree")

  result <- mcmc(model, z_init = c(0L, 0L, 1L),
                 psi_init = 0.5, n_iter = 100L,
                 psi_tune = 0.1, seed = 42L)

  eip <- result$edge_inclusion_probs(nug)
  expect_s3_class(eip, "data.frame")
  expect_named(eip, c("vertex1", "vertex2", "prob"))
  # Spanning tree on 3 vertices always has 2 edges, so probs should sum ~2
  expect_equal(nrow(eip), 3L)
  expect_true(all(eip$prob >= 0 & eip$prob <= 1))
})
