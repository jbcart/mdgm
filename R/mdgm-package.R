#' mdgm: Mixture of Directed Graphical Models
#'
#' Bayesian inference for discrete spatial random fields using mixtures of
#' directed graphical models (MDGM). The package provides MCMC sampling over
#' directed acyclic graph (DAG) structures---including spanning trees and
#' acyclic orientations---with optional emission models for hierarchical
#' observation processes.
#'
#' @section Model types:
#' \itemize{
#'   \item **Standalone**: The spatial field \eqn{z} is observed directly.
#'     The MCMC updates the DAG structure and dependence parameter \eqn{\psi}.
#'   \item **Hierarchical**: A latent spatial field \eqn{z} generates
#'     observations \eqn{y} through an emission distribution (e.g., Bernoulli).
#'     The MCMC additionally updates \eqn{z} and emission parameters.
#' }
#'
#' @section Key functions:
#' \itemize{
#'   \item Graph construction: [nug_from_edge_list()], [nug_from_adj_list()],
#'     [nug_from_adj_mat()]
#'   \item Model specification: `mdgm_model()`
#'   \item MCMC inference: `mcmc()`
#' }
#'
#' @references
#' Carter, J. B. and Calder, C. A. (2024). Mixture of Directed Graphical Models
#' for Discrete Spatial Random Fields. \doi{10.48550/arXiv.2406.15700}
#'
#' @docType package
#' @name mdgm
#' @keywords internal
#' @useDynLib mdgm, .registration = TRUE
"_PACKAGE"
