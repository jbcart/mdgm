# Create an undirected graph from an edge list

Create an undirected graph from an edge list

## Usage

``` r
nug_from_edge_list(n, edges, ...)
```

## Arguments

- n:

  Number of vertices.

- edges:

  Edge list as a matrix or data frame with two or three columns. The
  first two columns are vertex indices (1-indexed). Each undirected edge
  must appear twice (once per direction). An optional third column gives
  edge weights.

- ...:

  Additional arguments: `seed` (integer) for reproducible RNG.

## Value

A `NaturalUndirectedGraph` object.

## See also

Other graph-construction:
[`nug_from_adj_list()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_list.md),
[`nug_from_adj_mat()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_mat.md),
[`nug_from_grid()`](https://jbcart.github.io/mdgm/reference/nug_from_grid.md)

## Examples

``` r
# Triangle graph: vertices 1-2-3
edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
g <- nug_from_edge_list(3, edges)
g$nvertices()
#> [1] 3
g$neighbors(1)
#> [1] 2 3
```
