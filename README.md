
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mdgm

<!-- badges: start -->

<!-- badges: end -->

Full documentation: <https://jbcart.github.io/mdgm/>

**mdgm** (Mixture of Directed Graphical Models) provides Bayesian
inference for discrete spatial random fields. Instead of working with a
Markov random field and its intractable normalizing constant, the MDGM
defines a mixture over directed acyclic graphs (DAGs) compatible with an
undirected neighborhood graph. Each DAG admits a tractable likelihood,
avoiding the partition function entirely.

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
z <- c(1L, 1L, 0L, 0L,
       1L, 1L, 0L, 0L,
       0L, 0L, 1L, 1L,
       0L, 0L, 1L, 1L)
model <- mdgm_model(nug, dag_type = "spanning_tree")
result <- mcmc(model, z_init = z, psi_init = 0.5,
               n_iter = 2000L, psi_tune = 1.0, seed = 42L)
result$summary()
#> MDGM MCMC Results
#>   Vertices: 16, Colors: 2
#>   Iterations: 2000 (burnin: 0)
#>   Psi acceptance rate: 0.471
#>   Psi posterior mean: 0.8673 (sd: 0.5505)
#>   Psi R-hat: 1.0012, ESS: 259
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
