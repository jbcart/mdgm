# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

**mdgm** (Mixture of Directed Graphical Models) is an R package with a
C++20 core library. It implements graph data structures and sampling
algorithms for probabilistic graphical models. The C++ layer handles
computation; the R layer provides user-facing APIs via cpp11 bindings
and R6 classes. The paper is available at
<https://arxiv.org/html/2406.15700v4>.

## Build & Test Commands

### R Package

``` bash
R CMD build .
R CMD INSTALL .
```

### C++ Tests (GoogleTest)

``` bash
cmake -B build -DBUILD_TESTS=ON
cmake --build build
ctest --test-dir build
```

Run a single test suite:

``` bash
cmake --build build && ./build/tests/cpp/<test_binary>
```

### Regenerate cpp11 Bindings

After modifying `[[cpp11::register]]` annotated functions in
`src/R_*.cpp`:

``` r
cpp11::cpp_register()
```

This updates `src/cpp11.cpp` and `R/cpp11.R`.

## Architecture

### C++ Core (`inst/include/mdgm/` headers, `src/` implementations)

- **Graph Storage** (`graph_storage.hpp/cpp`) — Two sparse formats:
  `GraphCOO` (coordinate/edge list) and `GraphCSR` (compressed sparse
  row). COO is used for construction, CSR for efficient queries.
  Conversion and transposition supported.
- **NaturalUndirectedGraph** (`natural_undirected_graph.hpp/cpp`) —
  Undirected graph built on GraphCSR. Key algorithms:
  `SampleSpanningTree()` (Wilson’s, Aldous-Broder, or FastForward
  methods) and `SampleAcyclicOrientation()`. `GenerateRegularGraph()`
  creates 2D grid graphs.
- **DirectedAcyclicGraph** (`directed_acyclic_graph.hpp/cpp`) — DAG with
  both forward (CSR) and reverse (transposed) representations for
  parent/child queries.
- **RNG** (`rng.hpp`) — Wrapper around `std::mt19937_64` providing
  uniform, normal, discrete, and permutation sampling.

### R Interface (`R/`, `src/R_*.cpp`)

- **cpp11 bindings** (`src/R_natural_undirected_graph.cpp`,
  `src/R_utils.cpp`) — C++ objects are wrapped as R external pointers.
  Functions annotated with `[[cpp11::register]]` are auto-registered.
- **R6 class** (`R/natural_undirected_graph.R`) —
  `NaturalUndirectedGraph` R6 class wrapping the C++ object. User-facing
  constructors:
  [`nug_from_edge_list()`](https://jbcart.github.io/mdgm/reference/nug_from_edge_list.md),
  [`nug_from_adj_list()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_list.md),
  [`nug_from_adj_mat()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_mat.md).
- **Indexing convention**: R uses 1-indexed vectors; conversion to
  0-indexed C++ happens in the R binding layer
  (`R_natural_undirected_graph.cpp`).

### Tests (`tests/cpp/`)

C++ unit tests only (GoogleTest). Three test files mirror the core
classes: graph storage, undirected graphs, and DAGs.

## Key Conventions

- Directed edges follow the standard ordered pair: (parent, child) (the
  paper linked above currently has the reverse (child, parent))
- C++20 standard required (`CXX_STD = CXX20` in `src/Makevars`)
- Public C++ headers live in `inst/include/mdgm/` (R package convention
  for header-only/shared headers)
- Enum-based algorithm dispatch (e.g., `SpanningTreeMethod::kWilson`)
- Value semantics throughout C++ classes
