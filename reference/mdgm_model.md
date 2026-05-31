# Create an MDGM model (legacy interface)

**\[superseded\]**

`mdgm_model()` is superseded by
[`srf_model()`](https://jbcart.github.io/mdgm/reference/srf_model.md)
with [`mdgm()`](https://jbcart.github.io/mdgm/reference/mdgm.md)
configuration. It is kept for backwards compatibility.

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

An [SrfModel](https://jbcart.github.io/mdgm/reference/SrfModel.md)
object.
