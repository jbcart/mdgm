#' Sample from a Markov Random Field
#'
#' Draw a single sample from a binary (Ising) or multi-color (Potts) MRF on a
#' given graph.
#'
#' @param nug A `NaturalUndirectedGraph` object (external pointer).
#' @param psi Spatial dependence parameter (positive real).
#' @param n_colors Number of colors/categories (default 2, binary).
#' @param method Sampling method: `"auto"` (CFTP for k=2, Swendsen-Wang for
#'   k>2), `"swendsen_wang"`, `"gibbs"`, or `"cftp"` (k=2 only).
#' @param n_sweeps Number of sweeps (CFTP block size for binary, SW/Gibbs
#'   sweeps for multi-color). Default 200.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return Integer vector of length `nvertices(nug)` with values in
#'   `0, 1, ..., n_colors - 1`.
#'
#' @export
sample_mrf <- function(nug, psi, n_colors = 2L, method = "auto",
                       n_sweeps = 200L, seed = NULL) {
  stopifnot(is.numeric(psi), length(psi) == 1L, psi > 0)
  n_colors <- as.integer(n_colors)
  n_sweeps <- as.integer(n_sweeps)
  method <- match.arg(method, c("auto", "swendsen_wang", "gibbs", "cftp"))

  if (!is.null(seed)) {
    rng <- rng_create_seed_cpp(as.integer(seed))
  } else {
    rng <- rng_create_cpp()
  }

  nug_ptr <- if (inherits(nug, "NaturalUndirectedGraph")) nug$get_ptr() else nug
  sample_mrf_cpp(nug_ptr, psi, n_colors, n_sweeps, method, rng)
}
