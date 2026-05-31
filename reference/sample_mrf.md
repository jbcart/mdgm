# Sample from a Markov Random Field

Draw a single sample from a binary (Ising) or multi-color (Potts) MRF on
a given graph.

## Usage

``` r
sample_mrf(
  nug,
  psi,
  n_colors = 2L,
  method = "auto",
  n_sweeps = 200L,
  seed = NULL
)
```

## Arguments

- nug:

  A `NaturalUndirectedGraph` object (external pointer).

- psi:

  Spatial dependence parameter (positive real).

- n_colors:

  Number of colors/categories (default 2, binary).

- method:

  Sampling method: `"auto"` (CFTP for k=2, Swendsen-Wang for k\>2),
  `"swendsen_wang"`, `"gibbs"`, or `"cftp"` (k=2 only).

- n_sweeps:

  Number of sweeps (CFTP block size for binary, SW/Gibbs sweeps for
  multi-color). Default 200.

- seed:

  Optional integer seed for reproducibility.

## Value

Integer vector of length `nvertices(nug)` with values in
`0, 1, ..., n_colors - 1`.
