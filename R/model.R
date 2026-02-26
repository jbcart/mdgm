#' Spatial Random Field Model
#'
#' R6 class wrapping a C++ spatial random field model. Use [srf_model()] to
#' create instances.
#'
#' @importFrom R6 R6Class
#' @export
SrfModel <- R6::R6Class(
  classname = "SrfModel",
  cloneable = FALSE,
  public = list(
    #' @description Create a new SrfModel. Use [srf_model()] instead.
    #' @param ptr External pointer to a C++ Model object.
    #' @param model_type Character string: `"mdgm"` or `"mrf"`.
    initialize = function(ptr, model_type = "mdgm") {
      stopifnot(inherits(ptr, "externalptr"))
      private$.model <- ptr
      private$.model_type <- model_type
    },

    #' @description Check if the model has an emission layer.
    #' @return Logical.
    has_emission = function() {
      model_has_emission_cpp(private$.model)
    },

    #' @description Get the number of vertices.
    #' @return Integer.
    nvertices = function() {
      model_nvertices_cpp(private$.model)
    },

    #' @description Get the number of colors (categories).
    #' @return Integer.
    ncolors = function() {
      model_ncolors_cpp(private$.model)
    },

    #' @description Get the emission family type.
    #' @return Character string (`"bernoulli"`, `"gaussian"`, `"poisson"`) or
    #'   `NULL` for standalone models.
    emission_type = function() {
      model_emission_type_cpp(private$.model)
    },

    #' @description Get the spatial model type.
    #' @return Character string: `"mdgm"` or `"mrf"`.
    model_type = function() {
      private$.model_type
    },

    #' @description Get the internal C++ pointer. For internal use only.
    #' @return External pointer.
    #' @keywords internal
    get_ptr = function() {
      private$.model
    }
  ),
  private = list(
    .model = NULL,
    .model_type = NULL
  )
)

#' @rdname SrfModel
#' @usage NULL
#' @export
MdgmModel <- SrfModel

#' MDGM spatial configuration
#'
#' Creates a configuration object for a mixture of directed graphical models
#' spatial field. Used as the `spatial` argument to [srf_model()].
#'
#' @param dag_type DAG construction type: `"spanning_tree"`,
#'   `"acyclic_orientation"`, or `"rooted"`.
#' @return A list of class `mdgm_config`.
#' @examples
#' mdgm(dag_type = "spanning_tree")
#' @export
mdgm <- function(dag_type = c("spanning_tree",
                               "acyclic_orientation",
                               "rooted")) {
  dag_type <- match.arg(dag_type)
  structure(list(dag_type = dag_type), class = "mdgm_config")
}

#' MRF spatial configuration
#'
#' Creates a configuration object for a Markov random field (Potts model)
#' spatial field. Used as the `spatial` argument to [srf_model()].
#'
#' @param method Inference method for psi: `"exchange"` (exact, using the
#'   exchange algorithm) or `"pseudo_likelihood"` (approximate).
#' @param n_aux_sweeps Number of auxiliary Gibbs sweeps for the exchange
#'   algorithm (default 200). Ignored for pseudo-likelihood.
#' @return A list of class `mrf_config`.
#' @examples
#' mrf(method = "exchange")
#' mrf(method = "pseudo_likelihood")
#' @export
mrf <- function(method = c("exchange", "pseudo_likelihood"),
                n_aux_sweeps = 200L) {
  method <- match.arg(method)
  n_aux_sweeps <- as.integer(n_aux_sweeps)
  structure(list(method = method, n_aux_sweeps = n_aux_sweeps),
            class = "mrf_config")
}

#' Create a spatial random field model
#'
#' Constructs a spatial random field model for MCMC inference. The spatial
#' component can be either a mixture of directed graphical models (MDGM) or
#' a Markov random field (MRF), specified via the `spatial` argument using
#' the [mdgm()] or [mrf()] configuration helpers.
#'
#' @param nug A `NaturalUndirectedGraph` object defining the spatial structure.
#' @param spatial A spatial configuration object created by [mdgm()] or
#'   [mrf()].
#' @param emission Optional emission family for hierarchical models:
#'   `"bernoulli"`, `"gaussian"`, or `"poisson"`. If `NULL` (default),
#'   creates a standalone model where the spatial field is observed directly.
#' @param n_colors Number of categories for the spatial field (default 2).
#' @return An [SrfModel] object.
#' @examples
#' edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
#' nug <- nug_from_edge_list(3, edges, seed = 42L)
#'
#' # MDGM model
#' m1 <- srf_model(nug, spatial = mdgm(dag_type = "spanning_tree"))
#'
#' # MRF model with exchange algorithm
#' m2 <- srf_model(nug, spatial = mrf(method = "exchange"),
#'                 emission = "bernoulli")
#' @export
srf_model <- function(nug,
                      spatial = mdgm(),
                      emission = NULL,
                      n_colors = 2L) {
  stopifnot(inherits(nug, "NaturalUndirectedGraph"))
  n_colors <- as.integer(n_colors)
  nug_ptr <- nug$get_ptr()

  if (inherits(spatial, "mdgm_config")) {
    if (is.null(emission)) {
      ptr <- model_create_standalone_cpp(nug_ptr, spatial$dag_type, n_colors)
    } else {
      emission <- match.arg(emission, c("bernoulli", "gaussian", "poisson"))
      ptr <- model_create_hierarchical_cpp(
        nug_ptr, spatial$dag_type, n_colors, emission
      )
    }
    SrfModel$new(ptr, model_type = "mdgm")
  } else if (inherits(spatial, "mrf_config")) {
    if (is.null(emission)) {
      ptr <- model_create_mrf_standalone_cpp(
        nug_ptr, spatial$method, n_colors, spatial$n_aux_sweeps
      )
    } else {
      emission <- match.arg(emission, c("bernoulli", "gaussian", "poisson"))
      ptr <- model_create_mrf_hierarchical_cpp(
        nug_ptr, spatial$method, n_colors, emission, spatial$n_aux_sweeps
      )
    }
    SrfModel$new(ptr, model_type = "mrf")
  } else {
    stop("spatial must be created by mdgm() or mrf()")
  }
}

#' Create an MDGM model (legacy interface)
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#' `mdgm_model()` is superseded by [srf_model()] with [mdgm()] configuration.
#' It is kept for backwards compatibility.
#'
#' @param nug A `NaturalUndirectedGraph` object defining the spatial structure.
#' @param dag_type DAG construction type: `"spanning_tree"`,
#'   `"acyclic_orientation"`, or `"rooted"`.
#' @param n_colors Number of categories for the spatial field (default 2).
#' @param emission Optional emission family for hierarchical models:
#'   `"bernoulli"`, `"gaussian"`, or `"poisson"`. If `NULL` (default),
#'   creates a standalone model where the spatial field is observed directly.
#' @return An [SrfModel] object.
#' @export
mdgm_model <- function(nug,
                        dag_type = c("spanning_tree",
                                     "acyclic_orientation",
                                     "rooted"),
                        n_colors = 2L,
                        emission = NULL) {
  dag_type <- match.arg(dag_type)
  srf_model(nug, spatial = mdgm(dag_type = dag_type),
            emission = emission, n_colors = n_colors)
}
