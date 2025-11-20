#' @keywords internal
#' @importFrom R6 R6Class
UndirectedGraph <- R6::R6Class(
  classname = "UndirectedGraph",
  cloneable = FALSE,
  public = list(
    initialize = function(n, row_idx, col_idx, weights, ...) {
      args <- list(...)
      stopifnot(length(row_idx) == length(col_idx))
      if (any(row_idx < 1) || any(col_idx < 1) || any(row_idx > n) || any(col_idx > n)) {
        stop("Edge indices must be in the range [1, n]")
      }
      if (!is.integer(n)) {
        n <- as.integer(n)
      }
      if (!is.integer(row_idx)) {
        row <- as.integer(row_idx)
      }
      if (!is.integer(col_idx)) {
        col <- as.integer(col_idx)
      }
      private$.graph <- undirected_graph_create_cpp(n, row_idx, col_idx, weights)
      
      if ("rng" %in% names(args)) {
        stopifnot(inherits(args$rng, "externalptr"))
        private$.rng <- args$rng
      } else if ("seed" %in% names(args)) {
        private$.rng <- rng_create_seed_cpp(as.integer(args$seed))
      } else {
        private$.rng <- rng_create_cpp()
      }
    },
    
    nvertices = function() {
      undirected_graph_nvertices_cpp(private$.graph)
    },

    nedges = function() {
      undirected_graph_nedges_cpp(private$.graph)
    },

    neighbors = function(v) {
      undirected_graph_neighbors_cpp(private$.graph, as.integer(v))
    },

    sample_spanning_tree = function(method = c("wilson", "aldous_broder", "hybrid", "fast_forward"),
                                    k = 1000) {
      method <- match.arg(method)
      undirected_graph_sample_spanning_tree_cpp(private$.graph, method, as.integer(k), private$.rng)
    }
  ),
  private = list(
    .graph = NULL,
    .rng = NULL
  )
)

#' Creates an undirected graph from an edge list
#' 
#' @param n Number of vertices
#' @param edges Edge list as a matrix or data frame with two or three columns. The first two columns
#' represent the vertex indices of each edge. Edges must be represented twice, once for each
#' direction. The optional third column represents edge weights.
#' @param ... Additional arguments passed to the `UndirectedGraph` constructor for setting up the 
#' C++ mdgm::RNG. Can set the `seed` as an integer. 
#' @return An `UndirectedGraph` object.
#' @export
ug_from_edge_list <- function(n, edges, ...) {
  stopifnot(is.matrix(edges) || is.data.frame(edges))
  stopifnot(ncol(edges) == 2 || ncol(edges) == 3)
  row_idx <- as.integer(edges[, 1])
  col_idx <- as.integer(edges[, 2])
  if (ncol(edges) == 3) {
    weights <- as.numeric(edges[, 3])
  } else {
    weights <- rep(1.0, nrow(edges))
  }
  UndirectedGraph$new(n, row_idx, col_idx, weights, ...)
}

#' Creates an undirected graph from an adjacency list
#' 
#' @param adj_list Adjacency list as a list of integer vectors. The i-th element of the list
#' contains the neighbors of vertex i.
#' @return An `UndirectedGraph` object.
#' @param ... Additional arguments passed to the `UndirectedGraph` constructor for setting up the 
#' C++ mdgm::RNG. Can set the `seed` as an integer.
#' @export
ug_from_adj_list <- function(adj_list, ...) {
  n <- length(adj_list)
  edges <- Reduce(rbind, sapply(seq_along(adj_list), function(i) cbind(i, adj_list[[i]])))
  ug_from_edge_list(n, edges, ...)
}

#' Creates an undirected graph from an adjacency matrix
#'
#' @param adj_mat Adjacency matrix as a square symmetric matrix. Nonzero entries indicate the edges 
#' and edge weights. Self-loops are not allowed.
#' @param ... Additional arguments passed to the `UndirectedGraph` constructor for setting up the 
#' C++ mdgm::RNG. Can set the `seed` as an integer.
#' @return An `UndirectedGraph` object.
#' @export
ug_from_adj_mat <- function(adj_mat, ...) {
  stopifnot(is.matrix(adj_mat))
  stopifnot(isSymmetric(adj_mat))
  if (any(diag(adj_mat) != 0)) {
    stop("Self-loops are not allowed")
  }
  n <- nrow(adj_mat)
  edges <- which(adj_mat != 0, arr.ind = TRUE)
  edges <- cbind(edges, adj_mat[edges])
  ug_from_edge_list(n, edges, ...)
}


