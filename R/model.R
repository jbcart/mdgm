#' MDGM Model
#'
#' R6 class wrapping a C++ MDGM model. Use [mdgm_model()] to create instances.
#'
#' @importFrom R6 R6Class
#' @export
MdgmModel <- R6::R6Class(
  classname = "MdgmModel",
  cloneable = FALSE,
  public = list(
    #' @description Create a new MdgmModel. Use [mdgm_model()] instead.
    #' @param ptr External pointer to a C++ Model object.
    initialize = function(ptr) {
      stopifnot(inherits(ptr, "externalptr"))
      private$.model <- ptr
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
    }
  ),
  private = list(
    .model = NULL
  )
)

#' Create an MDGM model
#'
#' Constructs a mixture of directed graphical models for a discrete spatial
#' random field defined on an undirected graph.
#'
#' @param nug A `NaturalUndirectedGraph` object defining the spatial structure.
#' @param dag_type DAG construction type: `"spanning_tree"`,
#'   `"acyclic_orientation"`, or `"rooted"`.
#' @param n_colors Number of categories for the spatial field (default 2).
#' @param emission Optional emission family for hierarchical models:
#'   `"bernoulli"`. If `NULL` (default), creates a standalone model where
#'   the spatial field is observed directly.
#' @return An [MdgmModel] object.
#' @examples
#' # Standalone model on a triangle graph
#' edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
#' nug <- nug_from_edge_list(3, edges, seed = 42L)
#' model <- mdgm_model(nug, dag_type = "spanning_tree")
#' model$nvertices()
#' model$has_emission()
#' @export
mdgm_model <- function(nug,
                        dag_type = c("spanning_tree",
                                     "acyclic_orientation",
                                     "rooted"),
                        n_colors = 2L,
                        emission = NULL) {
  stopifnot(inherits(nug, "NaturalUndirectedGraph"))
  dag_type <- match.arg(dag_type)
  n_colors <- as.integer(n_colors)

  nug_ptr <- nug$.__enclos_env__$private$.graph

  if (is.null(emission)) {
    ptr <- model_create_standalone_cpp(nug_ptr, dag_type, n_colors)
  } else {
    emission <- match.arg(emission, c("bernoulli"))
    ptr <- model_create_hierarchical_cpp(
      nug_ptr, dag_type, n_colors, emission
    )
  }

  MdgmModel$new(ptr)
}
