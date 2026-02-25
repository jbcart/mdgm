# mdgm: Mixture of Directed Graphical Models

Bayesian inference for discrete spatial random fields using mixtures of
directed graphical models (MDGM). The package provides MCMC sampling
over directed acyclic graph (DAG) structures—including spanning trees
and acyclic orientations—with optional emission models for hierarchical
observation processes.

## Model types

- **Standalone**: The spatial field \\z\\ is observed directly. The MCMC
  updates the DAG structure and dependence parameter \\\psi\\.

- **Hierarchical**: A latent spatial field \\z\\ generates observations
  \\y\\ through an emission distribution (e.g., Bernoulli). The MCMC
  additionally updates \\z\\ and emission parameters.

## Key functions

- Graph construction:
  [`nug_from_edge_list()`](https://jbcart.github.io/mdgm/reference/nug_from_edge_list.md),
  [`nug_from_adj_list()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_list.md),
  [`nug_from_adj_mat()`](https://jbcart.github.io/mdgm/reference/nug_from_adj_mat.md)

- Model specification:
  [`mdgm_model()`](https://jbcart.github.io/mdgm/reference/mdgm_model.md)

- MCMC inference:
  [`mcmc()`](https://jbcart.github.io/mdgm/reference/mcmc.md)

## References

Carter, J. B. and Calder, C. A. (2024). Mixture of Directed Graphical
Models for Discrete Spatial Random Fields.
[doi:10.48550/arXiv.2406.15700](https://doi.org/10.48550/arXiv.2406.15700)

## See also

Useful links:

- <https://github.com/jbcart/mdgm>

- Report bugs at <https://github.com/jbcart/mdgm/issues>

## Author

**Maintainer**: Brandon Carter <brandon.carter@duke.edu>
