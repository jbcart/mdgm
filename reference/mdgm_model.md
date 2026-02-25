# Create an MDGM model

Constructs a mixture of directed graphical models for a discrete spatial
random field defined on an undirected graph.

## Usage

``` r
mdgm_model(
  nug,
  dag_type = c("spanning_tree", "acyclic_orientation", "rooted"),
  n_colors = 2L,
  emission = NULL
)
```

## Arguments

- nug:

  A `NaturalUndirectedGraph` object defining the spatial structure.

- dag_type:

  DAG construction type: `"spanning_tree"`, `"acyclic_orientation"`, or
  `"rooted"`.

- n_colors:

  Number of categories for the spatial field (default 2).

- emission:

  Optional emission family for hierarchical models: `"bernoulli"`,
  `"gaussian"`, or `"poisson"`. If `NULL` (default), creates a
  standalone model where the spatial field is observed directly.

## Value

An [MdgmModel](https://jbcart.github.io/mdgm/reference/MdgmModel.md)
object.

## Examples

``` r
# Standalone model on a triangle graph
edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
nug <- nug_from_edge_list(3, edges, seed = 42L)
model <- mdgm_model(nug, dag_type = "spanning_tree")
model$nvertices()
#> [1] 3
model$has_emission()
#> [1] FALSE
```
