#' Run MCMC inference for an MDGM model
#'
#' Performs Markov chain Monte Carlo sampling to estimate the posterior
#' distribution of the spatial field, DAG structure, dependence parameter,
#' and (for hierarchical models) emission parameters.
#'
#' @param model An [MdgmModel] object created by [mdgm_model()].
#' @param y Observation data. For hierarchical models, a list of integer
#'   vectors where `y[[i]]` contains the observations for vertex `i`.
#'   Vertices with no observations should have `integer(0)`. For standalone
#'   models, pass `NULL` (default).
#' @param z_init Initial color assignment, an integer vector of length `n`
#'   with values in `0, ..., n_colors - 1`.
#' @param psi_init Initial value for the dependence parameter (positive).
#' @param eta_init Initial emission parameters. For Bernoulli, a numeric
#'   vector of length `n_colors` (e.g., `c(0.2, 0.8)`). For Gaussian,
#'   a vector of length `2 * n_colors` with means then standard deviations
#'   (e.g., `c(mu_1, mu_2, sigma_1, sigma_2)`). For Poisson, a vector of
#'   length `n_colors` with rate parameters. Ignored for standalone models.
#' @param n_iter Number of MCMC iterations (default 1000).
#' @param psi_tune Standard deviation of the normal MH proposal for psi
#'   (default 0.1).
#' @param emission_prior_params Prior hyperparameters for emission
#'   parameters. For Bernoulli: `c(a, b)` for `Beta(a, b)` prior.
#'   For Gaussian: `c(mu_0, kappa_0, alpha_0, beta_0)` for
#'   Normal-InverseGamma prior. For Poisson: `c(alpha_0, beta_0)` for
#'   `Gamma(alpha_0, beta_0)` prior. Default `c(1, 1)`.
#' @param seed Optional integer seed for reproducibility.
#' @param nug Optional `NaturalUndirectedGraph` object. If provided, stored
#'   in the result for use by `edge_inclusion_probs()` and `plot()`.
#' @return An [MdgmResult] object containing posterior samples.
#' @seealso [mdgm_model()] for model construction, [MdgmResult] for
#'   accessing results.
#' @examples
#' \dontrun{
#' # Create a small grid graph
#' A <- matrix(0, 4, 4)
#' A[1, 2] <- A[2, 1] <- A[2, 3] <- A[3, 2] <- 1
#' A[3, 4] <- A[4, 3] <- 1
#' nug <- nug_from_adj_mat(A, seed = 42L)
#'
#' # Standalone model
#' model <- mdgm_model(nug, dag_type = "spanning_tree")
#' result <- mcmc(model, z_init = c(0L, 0L, 1L, 1L),
#'                psi_init = 0.5, n_iter = 100L)
#' result$summary()
#' }
#' @export
mcmc <- function(model, y = NULL, z_init, psi_init,
                 eta_init = numeric(0),
                 n_iter = 1000L, psi_tune = 0.1,
                 emission_prior_params = c(1, 1),
                 seed = NULL,
                 nug = NULL) {
  stopifnot(inherits(model, "MdgmModel"))

  n <- model$nvertices()
  z_init <- as.integer(z_init)
  stopifnot(length(z_init) == n)

  # Build flat observation storage
  if (is.null(y)) {
    obs_data <- integer(0)
    obs_ptr <- rep(0L, n + 1L)
  } else {
    stopifnot(is.list(y), length(y) == n)
    obs_data <- unlist(lapply(y, as.integer))
    lens <- vapply(y, length, integer(1))
    obs_ptr <- c(0L, cumsum(lens))
  }

  # Access internal model pointer
  model_ptr <- model$get_ptr()

  # Create RNG
  if (is.null(seed)) {
    rng_ptr <- rng_create_cpp()
  } else {
    rng_ptr <- rng_create_seed_cpp(as.integer(seed))
  }

  raw <- run_mcmc_cpp(
    model_ptr,
    as.integer(obs_data), as.integer(obs_ptr),
    z_init, as.double(psi_init),
    as.double(eta_init),
    as.integer(n_iter), as.double(psi_tune),
    as.double(emission_prior_params),
    rng_ptr
  )

  MdgmResult$new(raw,
                  emission_type = model$emission_type(),
                  nug = nug)
}
