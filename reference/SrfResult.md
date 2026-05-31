# MCMC Result

MCMC Result

MCMC Result

## Details

R6 class containing posterior samples from
[`mcmc()`](https://jbcart.github.io/mdgm/reference/mcmc.md). Access
samples via the accessor methods below.

## Methods

### Public methods

- [`SrfResult$new()`](#method-SrfResult-new)

- [`SrfResult$z()`](#method-SrfResult-z)

- [`SrfResult$alloc()`](#method-SrfResult-alloc)

- [`SrfResult$sufficient_stat()`](#method-SrfResult-sufficient_stat)

- [`SrfResult$z_map()`](#method-SrfResult-z_map)

- [`SrfResult$psi()`](#method-SrfResult-psi)

- [`SrfResult$emission_params()`](#method-SrfResult-emission_params)

- [`SrfResult$dag()`](#method-SrfResult-dag)

- [`SrfResult$edge_inclusion_probs()`](#method-SrfResult-edge_inclusion_probs)

- [`SrfResult$acceptance_rates()`](#method-SrfResult-acceptance_rates)

- [`SrfResult$diagnostics()`](#method-SrfResult-diagnostics)

- [`SrfResult$summary()`](#method-SrfResult-summary)

- [`SrfResult$plot()`](#method-SrfResult-plot)

------------------------------------------------------------------------

### Method `new()`

Create a new SrfResult. Use
[`mcmc()`](https://jbcart.github.io/mdgm/reference/mcmc.md) instead.

#### Usage

    SrfResult$new(raw, emission_type = NULL, nug = NULL, model_type = "mdgm")

#### Arguments

- `raw`:

  List of raw MCMC output from C++.

- `emission_type`:

  Character string or NULL.

- `nug`:

  A `NaturalUndirectedGraph` or NULL.

- `model_type`:

  Character string: `"mdgm"` or `"mrf"`.

------------------------------------------------------------------------

### Method `z()`

Get the latent field samples.

#### Usage

    SrfResult$z(iteration = NULL)

#### Arguments

- `iteration`:

  Optional integer iteration to extract (1-indexed). If `NULL`, returns
  the full `n x J` matrix.

#### Returns

Integer matrix (`n x J`) or integer vector (length `n`).

------------------------------------------------------------------------

### Method `alloc()`

Get allocation counts (always available).

#### Usage

    SrfResult$alloc()

#### Returns

Integer matrix (`n x k`) of per-site class counts across all iterations.

------------------------------------------------------------------------

### Method `sufficient_stat()`

Get per-iteration sufficient statistic T(z) (same-color edge count).

#### Usage

    SrfResult$sufficient_stat()

#### Returns

Numeric vector of length `J`.

------------------------------------------------------------------------

### Method `z_map()`

Get the joint MAP configuration.

#### Usage

    SrfResult$z_map()

#### Returns

Named list with `z` (integer vector), `log_posterior` (scalar), and
`iteration` (1-indexed integer).

------------------------------------------------------------------------

### Method `psi()`

Get the psi (dependence parameter) samples.

#### Usage

    SrfResult$psi()

#### Returns

Numeric vector of length `J`.

------------------------------------------------------------------------

### Method `emission_params()`

Get the emission parameters.

#### Usage

    SrfResult$emission_params(iteration = NULL)

#### Arguments

- `iteration`:

  Optional integer iteration to extract (1-indexed).

#### Returns

Named list appropriate to emission type, or `NULL` for standalone
models.

------------------------------------------------------------------------

### Method `dag()`

Get the DAG parent vector samples.

#### Usage

    SrfResult$dag(iteration = NULL)

#### Arguments

- `iteration`:

  Optional integer iteration to extract (1-indexed). If `NULL`, returns
  the full `n x J` matrix.

#### Returns

Integer matrix (`n x J`) or integer vector (length `n`), or `NULL` for
MRF models (which have no DAG structure). Values are 1-indexed parent
vertex IDs; `NA` indicates a root.

------------------------------------------------------------------------

### Method `edge_inclusion_probs()`

Compute edge inclusion probabilities.

#### Usage

    SrfResult$edge_inclusion_probs(nug = NULL, burnin = 0L)

#### Arguments

- `nug`:

  A `NaturalUndirectedGraph` object. If `NULL`, uses the graph stored at
  construction time.

- `burnin`:

  Number of initial iterations to discard (default 0).

#### Returns

A data frame with columns `vertex1`, `vertex2`, and `prob`.

------------------------------------------------------------------------

### Method `acceptance_rates()`

Get acceptance rates for MH steps.

#### Usage

    SrfResult$acceptance_rates()

#### Returns

A named numeric vector with `psi` and `graph` rates.

------------------------------------------------------------------------

### Method `diagnostics()`

Compute MCMC diagnostics.

#### Usage

    SrfResult$diagnostics(burnin = 0L)

#### Arguments

- `burnin`:

  Number of initial iterations to discard (default 0).

#### Returns

A named list with `rhat` and `ess` (effective sample size) for psi and
emission parameters.

------------------------------------------------------------------------

### Method [`summary()`](https://rdrr.io/r/base/summary.html)

Print a summary of the MCMC results.

#### Usage

    SrfResult$summary(burnin = 0L)

#### Arguments

- `burnin`:

  Number of initial iterations to discard (default 0).

#### Returns

Invisible named list with summary statistics.

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Plot MCMC trace plots and diagnostics.

#### Usage

    SrfResult$plot(burnin = 0L, which = "trace")

#### Arguments

- `burnin`:

  Number of initial iterations to discard (default 0).

- `which`:

  Character vector of what to plot: `"trace"` for trace plots of psi and
  emission params, `"edge_inclusion"` for edge inclusion probabilities
  (requires stored graph and igraph), `"posterior_field"` for posterior
  mode of the latent field (requires stored graph).

#### Returns

Invisible list of ggplot/igraph plot objects.
