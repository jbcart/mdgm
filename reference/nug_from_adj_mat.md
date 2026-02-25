# Create an undirected graph from an adjacency matrix

Create an undirected graph from an adjacency matrix

## Usage

``` r
nug_from_adj_mat(adj_mat, ...)
```

## Arguments

- adj_mat:

  A square symmetric matrix. Nonzero entries indicate edges; their
  values are used as edge weights. Self-loops are not allowed.

- ...:

  Additional arguments: `seed` (integer) for reproducible RNG.

## Value

A `NaturalUndirectedGraph` object.

## See also

Other graph-construction:
[`nug_from_adj_list()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_list.md),
[`nug_from_edge_list()`](https://jbcart.github.io/mdgm/reference/nug_from_edge_list.md),
[`nug_from_grid()`](https://jbcart.github.io/mdgm/reference/nug_from_grid.md)

## Examples

``` r
# 4-vertex path graph from adjacency matrix
A <- matrix(0, 4, 4)
A[1, 2] <- A[2, 1] <- 1
A[2, 3] <- A[3, 2] <- 1
A[3, 4] <- A[4, 3] <- 1
g <- nug_from_adj_mat(A)
g$nedges()
#> [1] 3
```
