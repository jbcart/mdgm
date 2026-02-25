# Create a regular 2D grid graph

Generates a grid graph with the specified dimensions and neighborhood
order. Order 1 gives rook adjacency (4-connected), order 2 gives queen
adjacency (8-connected).

## Usage

``` r
nug_from_grid(nrows, ncols, order = 1L, ...)
```

## Arguments

- nrows:

  Number of rows in the grid.

- ncols:

  Number of columns in the grid.

- order:

  Neighborhood order: `1` for rook (default) or `2` for queen.

- ...:

  Additional arguments: `seed` (integer) for reproducible RNG.

## Value

A `NaturalUndirectedGraph` object.

## See also

Other graph-construction:
[`nug_from_adj_list()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_list.md),
[`nug_from_adj_mat()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_mat.md),
[`nug_from_edge_list()`](https://jbcart.github.io/mdgm/reference/nug_from_edge_list.md)

## Examples

``` r
# 4x4 rook grid
g <- nug_from_grid(4, 4)
g$nvertices()
#> [1] 16
g$nedges()
#> [1] 24
```
