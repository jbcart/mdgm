# Create a spatial random field model

Constructs a spatial random field model for MCMC inference. The spatial
component can be either a mixture of directed graphical models (MDGM) or
a Markov random field (MRF), specified via the `spatial` argument using
the [`mdgm()`](https://jbcart.github.io/mdgm/reference/mdgm.md) or
[`mrf()`](https://jbcart.github.io/mdgm/reference/mrf.md) configuration
helpers.

## Usage

``` r
srf_model(nug, spatial = mdgm(), emission = NULL, n_colors = 2L)
```

## Arguments

- nug:

  A `NaturalUndirectedGraph` object defining the spatial structure.

- spatial:

  A spatial configuration object created by
  [`mdgm()`](https://jbcart.github.io/mdgm/reference/mdgm.md) or
  [`mrf()`](https://jbcart.github.io/mdgm/reference/mrf.md).

- emission:

  Optional emission family for hierarchical models: `"bernoulli"`,
  `"gaussian"`, or `"poisson"`. If `NULL` (default), creates a
  standalone model where the spatial field is observed directly.

- n_colors:

  Number of categories for the spatial field (default 2).

## Value

An [SrfModel](https://jbcart.github.io/mdgm/reference/SrfModel.md)
object.

## Examples

``` r
edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
nug <- nug_from_edge_list(3, edges, seed = 42L)

# MDGM model
m1 <- srf_model(nug, spatial = mdgm(dag_type = "spanning_tree"))

# MRF model with exchange algorithm
m2 <- srf_model(nug, spatial = mrf(method = "exchange"),
                emission = "bernoulli")
```
