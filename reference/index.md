# Package index

## Graph Construction

Create undirected graphs from various representations.

- [`nug_from_edge_list()`](https://jbcart.github.io/mdgm/reference/nug_from_edge_list.md)
  : Create an undirected graph from an edge list
- [`nug_from_adj_list()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_list.md)
  : Create an undirected graph from an adjacency list
- [`nug_from_adj_mat()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_mat.md)
  : Create an undirected graph from an adjacency matrix
- [`nug_from_grid()`](https://jbcart.github.io/mdgm/reference/nug_from_grid.md)
  : Create a regular 2D grid graph

## Model & Inference

Specify spatial random field models and run MCMC.

- [`srf_model()`](https://jbcart.github.io/mdgm/reference/srf_model.md)
  : Create a spatial random field model
- [`mdgm()`](https://jbcart.github.io/mdgm/reference/mdgm.md) : mdgm:
  Mixture of Directed Graphical Models
- [`mrf()`](https://jbcart.github.io/mdgm/reference/mrf.md) : MRF
  spatial configuration
- [`mcmc()`](https://jbcart.github.io/mdgm/reference/mcmc.md) : Run MCMC
  inference for a spatial random field model
- [`SrfModel`](https://jbcart.github.io/mdgm/reference/SrfModel.md)
  [`MdgmModel`](https://jbcart.github.io/mdgm/reference/SrfModel.md) :
  Spatial Random Field Model
- [`SrfResult`](https://jbcart.github.io/mdgm/reference/SrfResult.md) :
  MCMC Result
- [`sample_mrf()`](https://jbcart.github.io/mdgm/reference/sample_mrf.md)
  : Sample from a Markov Random Field

## Superseded

Functions that have been replaced by newer alternatives.

- [`mdgm_model()`](https://jbcart.github.io/mdgm/reference/mdgm_model.md)
  **\[superseded\]** : Create an MDGM model (legacy interface)
