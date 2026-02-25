# MCMC Result

MCMC Result

MCMC Result

## Details

R6 class containing posterior samples from
[`mcmc()`](https://jbcart.github.io/mdgm/reference/mcmc.md). Access
samples via the accessor methods below.

## Methods

### Public methods

- [`MdgmResult$new()`](#method-MdgmResult-new)

- [`MdgmResult$z()`](#method-MdgmResult-z)

- [`MdgmResult$psi()`](#method-MdgmResult-psi)

- [`MdgmResult$emission_params()`](#method-MdgmResult-emission_params)

- [`MdgmResult$dag()`](#method-MdgmResult-dag)

- [`MdgmResult$edge_inclusion_probs()`](#method-MdgmResult-edge_inclusion_probs)

- [`MdgmResult$acceptance_rates()`](#method-MdgmResult-acceptance_rates)

- [`MdgmResult$diagnostics()`](#method-MdgmResult-diagnostics)

- [`MdgmResult$summary()`](#method-MdgmResult-summary)

- [`MdgmResult$plot()`](#method-MdgmResult-plot)

------------------------------------------------------------------------

### Method `new()`

Create a new MdgmResult. Use
[`mcmc()`](https://jbcart.github.io/mdgm/reference/mcmc.md) instead.

#### Usage

    MdgmResult$new(raw, emission_type = NULL, nug = NULL)

#### Arguments

- `raw`:

  List of raw MCMC output from C++.

- `emission_type`:

  Character string or NULL.

- `nug`:

  A `NaturalUndirectedGraph` or NULL.

------------------------------------------------------------------------

### Method `z()`

Get the latent field samples.

#### Usage

    MdgmResult$z(iteration = NULL)

#### Arguments

- `iteration`:

  Optional integer iteration to extract (1-indexed). If `NULL`, returns
  the full `n x J` matrix.

#### Returns

Integer matrix (`n x J`) or integer vector (length `n`).

------------------------------------------------------------------------

### Method `psi()`

Get the psi (dependence parameter) samples.

#### Usage

    MdgmResult$psi()

#### Returns

Numeric vector of length `J`.

------------------------------------------------------------------------

### Method `emission_params()`

Get the emission parameters.

#### Usage

    MdgmResult$emission_params(iteration = NULL)

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

    MdgmResult$dag(iteration = NULL)

#### Arguments

- `iteration`:

  Optional integer iteration to extract (1-indexed). If `NULL`, returns
  the full `n x J` matrix.

#### Returns

Integer matrix (`n x J`) or integer vector (length `n`). Values are
1-indexed parent vertex IDs; `NA` indicates a root.

------------------------------------------------------------------------

### Method `edge_inclusion_probs()`

Compute edge inclusion probabilities.

#### Usage

    MdgmResult$edge_inclusion_probs(nug = NULL, burnin = 0L)

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

    MdgmResult$acceptance_rates()

#### Returns

A named numeric vector with `psi` and `graph` rates.

------------------------------------------------------------------------

### Method `diagnostics()`

Compute MCMC diagnostics.

#### Usage

    MdgmResult$diagnostics(burnin = 0L)

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

    MdgmResult$summary(burnin = 0L)

#### Arguments

- `burnin`:

  Number of initial iterations to discard (default 0).

#### Returns

Invisible named list with summary statistics.

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Plot MCMC trace plots and diagnostics.

#### Usage

    MdgmResult$plot(burnin = 0L, which = "trace")

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
