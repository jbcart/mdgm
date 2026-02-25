# Internal R6 class wrapping a C++ undirected graph

Internal R6 class wrapping a C++ undirected graph

Internal R6 class wrapping a C++ undirected graph

## Details

Users should create graphs via
[`nug_from_edge_list()`](https://jbcart.github.io/mdgm/reference/nug_from_edge_list.md),
[`nug_from_adj_list()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_list.md),
or
[`nug_from_adj_mat()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_mat.md)
rather than using this class directly.

## Methods

### Public methods

- [`NaturalUndirectedGraph$new()`](#method-NaturalUndirectedGraph-new)

- [`NaturalUndirectedGraph$nvertices()`](#method-NaturalUndirectedGraph-nvertices)

- [`NaturalUndirectedGraph$nedges()`](#method-NaturalUndirectedGraph-nedges)

- [`NaturalUndirectedGraph$neighbors()`](#method-NaturalUndirectedGraph-neighbors)

- [`NaturalUndirectedGraph$weights()`](#method-NaturalUndirectedGraph-weights)

- [`NaturalUndirectedGraph$sample_spanning_tree()`](#method-NaturalUndirectedGraph-sample_spanning_tree)

- [`NaturalUndirectedGraph$get_ptr()`](#method-NaturalUndirectedGraph-get_ptr)

------------------------------------------------------------------------

### Method `new()`

Create a new undirected graph.

#### Usage

    NaturalUndirectedGraph$new(
      n = NULL,
      row_idx = NULL,
      col_idx = NULL,
      weights = NULL,
      ...
    )

#### Arguments

- `n`:

  Number of vertices.

- `row_idx`:

  Integer vector of source vertex indices (1-indexed).

- `col_idx`:

  Integer vector of target vertex indices (1-indexed).

- `weights`:

  Numeric vector of edge weights.

- `...`:

  Additional arguments: `seed` (integer) or `rng` (external pointer).

------------------------------------------------------------------------

### Method `nvertices()`

Get the number of vertices.

#### Usage

    NaturalUndirectedGraph$nvertices()

#### Returns

Integer.

------------------------------------------------------------------------

### Method `nedges()`

Get the number of edges.

#### Usage

    NaturalUndirectedGraph$nedges()

#### Returns

Integer.

------------------------------------------------------------------------

### Method `neighbors()`

Get the 1-indexed neighbors of a vertex.

#### Usage

    NaturalUndirectedGraph$neighbors(v)

#### Arguments

- `v`:

  Vertex index (1-indexed).

#### Returns

Integer vector of neighbor indices.

------------------------------------------------------------------------

### Method [`weights()`](https://rdrr.io/r/stats/weights.html)

Get the edge weights to a vertex's neighbors.

#### Usage

    NaturalUndirectedGraph$weights(v)

#### Arguments

- `v`:

  Vertex index (1-indexed).

#### Returns

Numeric vector of edge weights (same order as `neighbors(v)`).

------------------------------------------------------------------------

### Method `sample_spanning_tree()`

Sample a uniform spanning tree.

#### Usage

    NaturalUndirectedGraph$sample_spanning_tree(
      method = c("wilson", "aldous_broder"),
      k = 1000
    )

#### Arguments

- `method`:

  Sampling algorithm: `"wilson"` or `"aldous_broder"`.

- `k`:

  Number of steps (unused, reserved for future methods).

#### Returns

A list representing the spanning tree as a parent vector.

------------------------------------------------------------------------

### Method `get_ptr()`

Get the internal C++ pointer. For internal use only.

#### Usage

    NaturalUndirectedGraph$get_ptr()

#### Returns

External pointer.
