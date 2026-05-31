# MRF spatial configuration

Creates a configuration object for a Markov random field (Potts model)
spatial field. Used as the `spatial` argument to
[`srf_model()`](https://jbcart.github.io/mdgm/reference/srf_model.md).

## Usage

``` r
mrf(method = c("exchange", "pseudo_likelihood"), n_aux_sweeps = 200L)
```

## Arguments

- method:

  Inference method for psi: `"exchange"` (exact, using the exchange
  algorithm) or `"pseudo_likelihood"` (approximate).

- n_aux_sweeps:

  Number of auxiliary Gibbs sweeps for the exchange algorithm (default
  200). Ignored for pseudo-likelihood.

## Value

A list of class `mrf_config`.

## Examples

``` r
mrf(method = "exchange")
#> $method
#> [1] "exchange"
#> 
#> $n_aux_sweeps
#> [1] 200
#> 
#> attr(,"class")
#> [1] "mrf_config"
mrf(method = "pseudo_likelihood")
#> $method
#> [1] "pseudo_likelihood"
#> 
#> $n_aux_sweeps
#> [1] 200
#> 
#> attr(,"class")
#> [1] "mrf_config"
```
