
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mdgm

<!-- badges: start -->

<!-- badges: end -->

Full documentation: <https://jbcart.github.io/mdgm/>

**mdgm** provides Bayesian inference for discrete spatial random fields.
It supports two spatial model families on undirected graphs:

- **Mixture of Directed Graphical Models (MDGM)** — avoids the
  intractable MRF partition function by defining a mixture over
  compatible DAGs.
- **Markov Random Field (MRF)** — classical Potts/Ising model with
  inference via the exchange algorithm (exact) or pseudo-likelihood
  (approximate).

Both can be combined with emission distributions (Bernoulli, Gaussian,
Poisson) for hierarchical models where the spatial field is latent.

See [Carter and Calder (2024)](https://arxiv.org/abs/2406.15700) for the
full methodological details.

## Installation

You can install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("jbcart/mdgm")
```

A C++20 compiler is required.

## Quick start

``` r
library(mdgm)

# Build a 4x4 grid graph with rook adjacency
nug <- nug_from_grid(4, 4, seed = 42L)
nug$nvertices()
#> [1] 16
nug$nedges()
#> [1] 24

# Fit a standalone spanning-tree MDGM
z <- c(0L, 0L, 0L, 1L,
       0L, 0L, 1L, 1L,
       1L, 1L, 1L, 0L,
       1L, 1L, 0L, 0L)
model <- srf_model(nug, spatial = mdgm(dag_type = "spanning_tree"))
result <- mcmc(model, z_init = z, psi_init = 0.5,
               n_iter = 2000L, psi_tune = 1.0, seed = 42L)
result$summary()
#> MDGM MCMC Results
#>   Vertices: 16, Colors: 2
#>   Iterations: 2000 (burnin: 0)
#>   Psi acceptance rate: 0.462
#>   Psi posterior mean: 0.8447 (sd: 0.5658)
#>   Diagnostics:
#>     psi — R-hat: 1.0035, ESS: 308
```

### MRF example

``` r
# Fit an MRF with pseudo-likelihood
model_mrf <- srf_model(nug, spatial = mrf(method = "pseudo_likelihood"))
result_mrf <- mcmc(model_mrf, z_init = z, psi_init = 0.5,
                   n_iter = 2000L, psi_tune = 0.5, seed = 42L)
result_mrf$summary()
#> MRF MCMC Results
#>   Vertices: 16, Colors: 2
#>   Iterations: 2000 (burnin: 0)
#>   Psi acceptance rate: 0.699
#>   Psi posterior mean: 1.0303 (sd: 0.5364)
#>   Diagnostics:
#>     psi — R-hat: 1.0044, ESS: 110
```

See `vignette("mdgm")` for a full walkthrough including visualization
and posterior diagnostics.

## C++ unit tests

The C++ core ships with GoogleTest unit tests under `tests/cpp/`. Build
and run them with:

``` bash
cmake -B build -DBUILD_TESTS=ON
cmake --build build
ctest --test-dir build
```
