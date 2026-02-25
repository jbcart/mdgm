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
                 theta_init = c(0.3, 0.7),
                 n_iter = 50L,
                 psi_tune = 0.1,
                 emission_prior_params = c(1, 1),
                 seed = 123L)

  expect_s3_class(result, "MdgmResult")

  # Check emission_params
  ep <- result$emission_params()
  expect_type(ep, "list")
  expect_named(ep, "p")
  expect_equal(dim(ep$p), c(2L, 50L))

  # p values should be between 0 and 1
  expect_true(all(ep$p >= 0 & ep$p <= 1))
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

test_that("summary prints and returns structured data", {
  edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree")

  result <- mcmc(model, z_init = c(0L, 0L, 1L),
                 psi_init = 0.5, n_iter = 50L,
                 psi_tune = 0.1, seed = 42L)

  out <- expect_output(result$summary(), "MDGM MCMC Results")
  expect_type(out, "list")
  expect_true("psi_mean" %in% names(out))
  expect_true("diagnostics" %in% names(out))
})

test_that("diagnostics returns R-hat and ESS", {
  edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree",
                       emission = "bernoulli")
  y <- list(c(1L, 0L), c(0L, 0L), c(1L, 1L))

  result <- mcmc(model, y = y, z_init = c(0L, 0L, 1L),
                 psi_init = 0.5, theta_init = c(0.3, 0.7),
                 n_iter = 200L, psi_tune = 0.1,
                 emission_prior_params = c(1, 1), seed = 42L)

  diag <- result$diagnostics(burnin = 50L)
  expect_type(diag, "list")
  expect_true("psi" %in% names(diag))
  expect_true(is.finite(diag$psi$rhat))
  expect_true(diag$psi$ess > 0)
  expect_true("p" %in% names(diag))
  expect_length(diag$p$rhat, 2L)
})

test_that("standalone MCMC works with n_colors = 3", {
  edges <- rbind(
    c(1, 2), c(2, 1), c(2, 3), c(3, 2),
    c(1, 3), c(3, 1), c(3, 4), c(4, 3)
  )
  nug <- nug_from_edge_list(4, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree", n_colors = 3L)

  result <- mcmc(model,
                 z_init = c(0L, 1L, 2L, 0L),
                 psi_init = 0.5,
                 n_iter = 50L,
                 psi_tune = 0.1,
                 seed = 123L)

  z_mat <- result$z()
  expect_equal(dim(z_mat), c(4L, 50L))
  expect_true(all(z_mat %in% c(0L, 1L, 2L)))
})

test_that("hierarchical MCMC works with n_colors = 3 Bernoulli", {
  edges <- rbind(
    c(1, 2), c(2, 1), c(2, 3), c(3, 2),
    c(1, 3), c(3, 1), c(3, 4), c(4, 3)
  )
  nug <- nug_from_edge_list(4, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree",
                       n_colors = 3L, emission = "bernoulli")

  y <- list(c(1L, 0L, 1L), c(0L, 0L), c(1L, 1L, 1L), c(0L, 1L))

  result <- mcmc(model, y = y,
                 z_init = c(0L, 0L, 2L, 1L),
                 psi_init = 0.5,
                 theta_init = c(0.2, 0.5, 0.8),
                 n_iter = 50L,
                 psi_tune = 0.1,
                 emission_prior_params = c(1, 1),
                 seed = 123L)

  ep <- result$emission_params()
  p_mat <- ep$p
  expect_equal(dim(p_mat), c(3L, 50L))

  # All p in [0, 1]
  expect_true(all(p_mat >= 0 & p_mat <= 1))

  # Ordering: p[1,] < p[2,] < p[3,] for each iteration
  for (j in seq_len(ncol(p_mat))) {
    expect_true(p_mat[1, j] < p_mat[2, j],
                info = paste("iter", j))
    expect_true(p_mat[2, j] < p_mat[3, j],
                info = paste("iter", j))
  }
})

test_that("hierarchical MCMC works with Gaussian emission", {
  edges <- rbind(
    c(1, 2), c(2, 1), c(2, 3), c(3, 2),
    c(1, 3), c(3, 1), c(3, 4), c(4, 3)
  )
  nug <- nug_from_edge_list(4, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree",
                       emission = "gaussian")

  # Gaussian observations (encoded as integers)
  y <- list(c(1L, 2L, 3L), c(2L, 1L), c(8L, 9L, 10L), c(7L, 8L))

  result <- mcmc(model, y = y,
                 z_init = c(0L, 0L, 1L, 1L),
                 psi_init = 0.5,
                 theta_init = c(2.0, 8.0, 1.0, 1.0),
                 n_iter = 50L,
                 psi_tune = 0.1,
                 emission_prior_params = c(0, 10000, 2, 1),
                 seed = 123L)

  ep <- result$emission_params()
  expect_type(ep, "list")
  expect_true("mu" %in% names(ep))
  expect_true("sigma2" %in% names(ep))
})

test_that("hierarchical MCMC works with Poisson emission", {
  edges <- rbind(
    c(1, 2), c(2, 1), c(2, 3), c(3, 2),
    c(1, 3), c(3, 1), c(3, 4), c(4, 3)
  )
  nug <- nug_from_edge_list(4, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree",
                       emission = "poisson")

  y <- list(c(0L, 1L, 0L), c(1L, 0L), c(5L, 6L, 4L), c(3L, 5L))

  result <- mcmc(model, y = y,
                 z_init = c(0L, 0L, 1L, 1L),
                 psi_init = 0.5,
                 theta_init = c(1.0, 5.0),
                 n_iter = 50L,
                 psi_tune = 0.1,
                 emission_prior_params = c(1, 0.1),
                 seed = 123L)

  ep <- result$emission_params()
  expect_type(ep, "list")
  expect_named(ep, "lambda")
  # Lambda should be positive
  expect_true(all(ep$lambda > 0))
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

test_that("edge inclusion probabilities work with burnin", {
  edges <- rbind(
    c(1, 2), c(2, 1), c(2, 3), c(3, 2),
    c(1, 3), c(3, 1)
  )
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree")

  result <- mcmc(model, z_init = c(0L, 0L, 1L),
                 psi_init = 0.5, n_iter = 100L,
                 psi_tune = 0.1, seed = 42L)

  eip <- result$edge_inclusion_probs(nug, burnin = 50L)
  expect_equal(nrow(eip), 3L)
  expect_true(all(eip$prob >= 0 & eip$prob <= 1))

  # Burnin exceeding n_iter should error
  expect_error(result$edge_inclusion_probs(nug, burnin = 100L))
})

test_that("edge_inclusion_probs uses stored graph", {
  edges <- rbind(
    c(1, 2), c(2, 1), c(2, 3), c(3, 2),
    c(1, 3), c(3, 1)
  )
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree")

  result <- mcmc(model, z_init = c(0L, 0L, 1L),
                 psi_init = 0.5, n_iter = 50L,
                 psi_tune = 0.1, seed = 42L,
                 nug = nug)

  # Should work without passing nug

  eip <- result$edge_inclusion_probs()
  expect_s3_class(eip, "data.frame")
  expect_equal(nrow(eip), 3L)
})

test_that("standalone emission_params returns NULL", {
  edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2))
  nug <- nug_from_edge_list(3, edges, seed = 42L)
  model <- mdgm_model(nug, dag_type = "spanning_tree")

  result <- mcmc(model, z_init = c(0L, 0L, 1L),
                 psi_init = 0.5, n_iter = 20L,
                 psi_tune = 0.1, seed = 42L)

  expect_null(result$emission_params())
})
