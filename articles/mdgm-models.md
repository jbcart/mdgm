# MDGM Model Specification

## Background

Classical spatial models for discrete data typically use a Markov random
field (MRF), which specifies a joint distribution over an undirected
graph. However, the MRF likelihood involves an intractable normalizing
constant (the partition function), making exact Bayesian inference
computationally burdensome.

The **Mixture of Directed Graphical Models (MDGM)** is an alternative
that takes the same undirected graph as input but defines a mixture over
compatible directed acyclic graphs (DAGs). Each DAG admits a tractable
factorization, so the MDGM avoids the partition function entirely. The
joint distribution marginalizes over the DAG space:
$p(\mathbf{z} \mid {\mathbf{ξ}}) = \sum_{D}p(D)\, p(\mathbf{z} \mid D,{\mathbf{ξ}})$.

For full details, see [Carter and Calder
(2024)](https://arxiv.org/abs/2406.15700).

## The MDGM prior

Let $G = (V,E)$ be an undirected graph (the “natural undirected graph”
or NUG) encoding potential neighbor relationships. The MDGM places a
prior over DAGs $D$ that are compatible with $G$: every directed edge
$(u,v)$ in $D$ corresponds to an undirected edge $\{ u,v\}$ in $G$.

Three DAG constructions are supported:

### Spanning trees

A spanning tree $T$ of $G$ is a connected, acyclic subgraph containing
all vertices. Edges are directed from child to parent. The posterior
sampling uses Wilson’s algorithm with data-dependent edge weights:

$$w(u,v) = \exp(\psi \cdot \mathbf{1}\left\lbrack z_{u} = z_{v} \right\rbrack)$$

where $\psi > 0$ is the spatial dependence parameter and $z$ is the
color assignment. This provides a direct (non-MH) posterior sample of
the spanning tree.

### Acyclic orientations

An acyclic orientation assigns a direction to every edge in $G$ such
that no directed cycle exists. Equivalently, this is defined by a vertex
permutation $\sigma$: edge $\{ u,v\}$ is directed as $(u,v)$ if
$\sigma(u) > \sigma(v)$. The MCMC proposes new permutations and accepts
via a Metropolis-Hastings step based on the exact DAG log-likelihood
ratio.

### Rooted DAGs

A rooted DAG is constructed from $G$ by choosing a root vertex and
orienting edges via a breadth-first or depth-first traversal. The MCMC
updates the root via random walk proposals on $G$.

## Spatial field model

Given a DAG $D$, each vertex $i$ has a set of parents
$\text{pa}_{D}(i)$. The conditional distribution of the color $z_{i}$
given its parents is:

\$\$p(z_i = k \mid z\_{\text{pa}(i)}, \psi) \propto \exp\Bigl(\psi
\sum\_{j \in \text{pa}(i)} \mathbf{1}\[z_j = k\] + \alpha_k\Bigr)\$\$

where $\alpha_{k}$ are marginal log-probabilities (currently fixed at 0
for a uniform marginal).

## Standalone vs. hierarchical models

### Standalone model

In the standalone model, the spatial field $z$ is observed directly. The
MCMC updates only the DAG structure $D$ and the dependence parameter
$\psi$. This is useful when the data are categorical labels on a spatial
domain.

### Hierarchical model

In the hierarchical model, $z$ is a latent field and observations
$y_{i}$ are generated through an **emission distribution**:

$$y_{ij} \mid z_{i},\eta \sim f\left( y_{ij} \mid \eta_{z_{i}} \right)$$

Currently supported emission families:

- **Bernoulli**:
  $y_{ij} \mid z_{i} = k \sim \text{Bernoulli}\left( \eta_{k} \right)$,
  with identifiability constraint $\eta_{0} < \eta_{1} < \cdots$
  enforced via truncated Beta posterior sampling.
- **Gaussian**:
  $y_{ij} \mid z_{i} = k \sim \mathcal{N}\left( \mu_{k},\sigma_{k}^{2} \right)$,
  with identifiability constraint $\mu_{0} < \mu_{1} < \cdots$.
  Conjugate Normal-InverseGamma updates.
- **Poisson**:
  $y_{ij} \mid z_{i} = k \sim \text{Poisson}\left( \lambda_{k} \right)$,
  with identifiability constraint $\lambda_{0} < \lambda_{1} < \cdots$.
  Conjugate truncated Gamma updates.

The MCMC additionally updates $z$ (Gibbs scan over vertices) and the
emission parameters (conjugate updates). See the [Emission
Models](https://jbcart.github.io/mdgm/articles/emission-models.md)
vignette for worked examples.

## Prior specification

### Dependence parameter

The spatial dependence parameter $\psi > 0$ has a half-Cauchy prior:

$$p(\psi) = \frac{2}{\pi\left( 1 + \psi^{2} \right)},\quad\psi > 0$$

Updates use a Metropolis-Hastings random walk with a normal proposal.

### Emission parameters

- **Bernoulli**: Each $\eta_{k}$ has a $\text{Beta}(a,b)$ prior.
  `emission_prior_params = c(a, b)`.
- **Gaussian**: Each $\left( \mu_{k},\sigma_{k}^{2} \right)$ has a
  Normal-InverseGamma prior with hyperparameters
  $\left( \mu_{0},\kappa_{0},\alpha_{0},\beta_{0} \right)$.
  `emission_prior_params = c(mu_0, kappa_0, alpha_0, beta_0)`.
- **Poisson**: Each $\lambda_{k}$ has a
  $\text{Gamma}\left( \alpha_{0},\beta_{0} \right)$ prior (rate
  parameterization). `emission_prior_params = c(alpha_0, beta_0)`.

All conjugate posteriors use truncated sampling to enforce parameter
ordering for identifiability.

## MCMC algorithm

Each iteration of the MCMC sampler performs:

1.  **Update DAG** $D$ — Sample a new spanning tree (direct posterior
    sample via Wilson’s) or propose a new acyclic orientation/root (MH
    step).
2.  **Update** $z$ (hierarchical only) — Gibbs scan over vertices in
    random order. The full conditional combines the spatial prior with
    the emission likelihood.
3.  **Update** $\psi$ — Metropolis-Hastings with normal random walk
    proposal.
4.  **Update** $\eta$ (hierarchical only) — Conjugate posterior sampling
    with identifiability constraints.
