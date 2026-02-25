# Run MCMC inference for an MDGM model

Performs Markov chain Monte Carlo sampling to estimate the posterior
distribution of the spatial field, DAG structure, dependence parameter,
and (for hierarchical models) emission parameters.

## Usage

``` r
mcmc(
  model,
  y = NULL,
  z_init,
  psi_init,
  eta_init = numeric(0),
  n_iter = 1000L,
  psi_tune = 0.1,
  emission_prior_params = c(1, 1),
  seed = NULL,
  nug = NULL
)
```

## Arguments

- model:

  An [MdgmModel](https://jbcart.github.io/mdgm/reference/MdgmModel.md)
  object created by
  [`mdgm_model()`](https://jbcart.github.io/mdgm/reference/mdgm_model.md).

- y:

  Observation data. For hierarchical models, a list of integer vectors
  where `y[[i]]` contains the observations for vertex `i`. Vertices with
  no observations should have `integer(0)`. For standalone models, pass
  `NULL` (default).

- z_init:

  Initial color assignment, an integer vector of length `n` with values
  in `0, ..., n_colors - 1`.

- psi_init:

  Initial value for the dependence parameter (positive).

- eta_init:

  Initial emission parameters. For Bernoulli, a numeric vector of length
  `n_colors` (e.g., `c(0.2, 0.8)`). For Gaussian, a vector of length
  `2 * n_colors` with means then standard deviations (e.g.,
  `c(mu_1, mu_2, sigma_1, sigma_2)`). For Poisson, a vector of length
  `n_colors` with rate parameters. Ignored for standalone models.

- n_iter:

  Number of MCMC iterations (default 1000).

- psi_tune:

  Standard deviation of the normal MH proposal for psi (default 0.1).

- emission_prior_params:

  Prior hyperparameters for emission parameters. For Bernoulli:
  `c(a, b)` for `Beta(a, b)` prior. For Gaussian:
  `c(mu_0, kappa_0, alpha_0, beta_0)` for Normal-InverseGamma prior. For
  Poisson: `c(alpha_0, beta_0)` for `Gamma(alpha_0, beta_0)` prior.
  Default `c(1, 1)`.

- seed:

  Optional integer seed for reproducibility.

- nug:

  Optional `NaturalUndirectedGraph` object. If provided, stored in the
  result for use by `edge_inclusion_probs()` and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

An [MdgmResult](https://jbcart.github.io/mdgm/reference/MdgmResult.md)
object containing posterior samples.

## See also

[`mdgm_model()`](https://jbcart.github.io/mdgm/reference/mdgm_model.md)
for model construction,
[MdgmResult](https://jbcart.github.io/mdgm/reference/MdgmResult.md) for
accessing results.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a small grid graph
A <- matrix(0, 4, 4)
A[1, 2] <- A[2, 1] <- A[2, 3] <- A[3, 2] <- 1
A[3, 4] <- A[4, 3] <- 1
nug <- nug_from_adj_mat(A, seed = 42L)

# Standalone model
model <- mdgm_model(nug, dag_type = "spanning_tree")
result <- mcmc(model, z_init = c(0L, 0L, 1L, 1L),
               psi_init = 0.5, n_iter = 100L)
result$summary()
} # }
```
