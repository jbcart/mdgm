#' Sample from a Markov Random Field
#'
#' Draw a single sample from a binary (Ising) or multi-color (Potts) MRF on a
#' given graph. For binary fields (`n_colors = 2`), uses CFTP (Propp-Wilson)
#' for exact sampling; otherwise uses Gibbs sampling.
#'
#' @param nug A `NaturalUndirectedGraph` object (external pointer).
#' @param psi Spatial dependence parameter (positive real).
#' @param n_colors Number of colors/categories (default 2, binary).
#' @param n_sweeps Number of Gibbs sweeps (used as CFTP block size for binary,
#'   or total sweeps for Potts). Default 200.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return Integer vector of length `nvertices(nug)` with values in
#'   `0, 1, ..., n_colors - 1`.
#'
#' @export
sample_mrf <- function(nug, psi, n_colors = 2L, n_sweeps = 200L, seed = NULL) {
  stopifnot(is.numeric(psi), length(psi) == 1L, psi > 0)
  n_colors <- as.integer(n_colors)
  n_sweeps <- as.integer(n_sweeps)

  if (!is.null(seed)) {
    rng <- rng_create_seed_cpp(as.integer(seed))
  } else {
    rng <- rng_create_cpp()
  }

  nug_ptr <- if (inherits(nug, "NaturalUndirectedGraph")) nug$get_ptr() else nug
  sample_mrf_cpp(nug_ptr, psi, n_colors, n_sweeps, rng)
}
