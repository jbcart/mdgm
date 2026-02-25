# Create an undirected graph from an adjacency list

Create an undirected graph from an adjacency list

## Usage

``` r
nug_from_adj_list(adj_list, ...)
```

## Arguments

- adj_list:

  A list of integer vectors. The i-th element contains the 1-indexed
  neighbors of vertex i.

- ...:

  Additional arguments: `seed` (integer) for reproducible RNG.

## Value

A `NaturalUndirectedGraph` object.

## See also

Other graph-construction:
[`nug_from_adj_mat()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_mat.md),
[`nug_from_edge_list()`](https://jbcart.github.io/mdgm/reference/nug_from_edge_list.md),
[`nug_from_grid()`](https://jbcart.github.io/mdgm/reference/nug_from_grid.md)

## Examples

``` r
# Path graph: 1 -- 2 -- 3
adj <- list(c(2L), c(1L, 3L), c(2L))
g <- nug_from_adj_list(adj)
g$nvertices()
#> [1] 3
```
