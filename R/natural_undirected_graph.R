#' Internal R6 class wrapping a C++ undirected graph
#'
#' Users should create graphs via [nug_from_edge_list()], [nug_from_adj_list()],
#' or [nug_from_adj_mat()] rather than using this class directly.
#'
#' @keywords internal
#' @importFrom R6 R6Class
NaturalUndirectedGraph <- R6::R6Class(
  classname = "NaturalUndirectedGraph",
  cloneable = FALSE,
  public = list(
    #' @description Create a new undirected graph.
    #' @param n Number of vertices.
    #' @param row_idx Integer vector of source vertex indices (1-indexed).
    #' @param col_idx Integer vector of target vertex indices (1-indexed).
    #' @param weights Numeric vector of edge weights.
    #' @param ... Additional arguments: `seed` (integer) or `rng` (external pointer).
    initialize = function(n = NULL, row_idx = NULL, col_idx = NULL, weights = NULL, ...) {
      args <- list(...)
      if (".graph_ptr" %in% names(args)) {
        private$.graph <- args[[".graph_ptr"]]
      } else {
        stopifnot(!is.null(n), !is.null(row_idx), !is.null(col_idx), !is.null(weights))
        stopifnot(length(row_idx) == length(col_idx))
        if (any(row_idx < 1) || any(col_idx < 1) || any(row_idx > n) || any(col_idx > n)) {
          stop("Edge indices must be in the range [1, n]")
        }
        if (!is.integer(n)) {
          n <- as.integer(n)
        }
        if (!is.integer(row_idx)) {
          row_idx <- as.integer(row_idx)
        }
        if (!is.integer(col_idx)) {
          col_idx <- as.integer(col_idx)
        }
        private$.graph <- nug_create_cpp(n, row_idx, col_idx, weights)
      }
      if ("rng" %in% names(args)) {
        stopifnot(inherits(args$rng, "externalptr"))
        private$.rng <- args$rng
      } else if ("seed" %in% names(args)) {
        private$.rng <- rng_create_seed_cpp(as.integer(args$seed))
      } else {
        private$.rng <- rng_create_cpp()
      }
    },
    
    #' @description Get the number of vertices.
    #' @return Integer.
    nvertices = function() {
      nug_nvertices_cpp(private$.graph)
    },

    #' @description Get the number of edges.
    #' @return Integer.
    nedges = function() {
      nug_nedges_cpp(private$.graph)
    },

    #' @description Get the 1-indexed neighbors of a vertex.
    #' @param v Vertex index (1-indexed).
    #' @return Integer vector of neighbor indices.
    neighbors = function(v) {
      nug_neighbors_cpp(private$.graph, as.integer(v))
    },

    #' @description Get the edge weights to a vertex's neighbors.
    #' @param v Vertex index (1-indexed).
    #' @return Numeric vector of edge weights (same order as `neighbors(v)`).
    weights = function(v) {
      nug_neighbor_weights_cpp(private$.graph, as.integer(v))
    },

    #' @description Sample a uniform spanning tree.
    #' @param method Sampling algorithm: `"wilson"` or `"aldous_broder"`.
    #' @param k Number of steps (unused, reserved for future methods).
    #' @return A list representing the spanning tree as a parent vector.
    sample_spanning_tree = function(method = c("wilson", "aldous_broder"),
                                    k = 1000) {
      method <- match.arg(method)
      nug_sample_spanning_tree_cpp(private$.graph, method, as.integer(k), private$.rng)
    },

    #' @description Get the internal C++ pointer. For internal use only.
    #' @return External pointer.
    #' @keywords internal
    get_ptr = function() {
      private$.graph
    }
  ),
  private = list(
    .graph = NULL,
    .rng = NULL
  )
)

#' Create an undirected graph from an edge list
#'
#' @param n Number of vertices.
#' @param edges Edge list as a matrix or data frame with two or three columns.
#'   The first two columns are vertex indices (1-indexed). Each undirected edge
#'   must appear twice (once per direction). An optional third column gives edge
#'   weights.
#' @param ... Additional arguments: `seed` (integer) for reproducible RNG.
#' @return A `NaturalUndirectedGraph` object.
#' @family graph-construction
#' @examples
#' # Triangle graph: vertices 1-2-3
#' edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
#' g <- nug_from_edge_list(3, edges)
#' g$nvertices()
#' g$neighbors(1)
#' @export
nug_from_edge_list <- function(n, edges, ...) {
  stopifnot(is.matrix(edges) || is.data.frame(edges))
  stopifnot(ncol(edges) == 2 || ncol(edges) == 3)
  row_idx <- as.integer(edges[, 1])
  col_idx <- as.integer(edges[, 2])
  if (ncol(edges) == 3) {
    weights <- as.numeric(edges[, 3])
  } else {
    weights <- rep(1.0, nrow(edges))
  }
  NaturalUndirectedGraph$new(n, row_idx, col_idx, weights, ...)
}

#' Create an undirected graph from an adjacency list
#'
#' @param adj_list A list of integer vectors. The i-th element contains the
#'   1-indexed neighbors of vertex i.
#' @param ... Additional arguments: `seed` (integer) for reproducible RNG.
#' @return A `NaturalUndirectedGraph` object.
#' @family graph-construction
#' @examples
#' # Path graph: 1 -- 2 -- 3
#' adj <- list(c(2L), c(1L, 3L), c(2L))
#' g <- nug_from_adj_list(adj)
#' g$nvertices()
#' @export
nug_from_adj_list <- function(adj_list, ...) {
  n <- length(adj_list)
  edges <- do.call(rbind, lapply(seq_along(adj_list), function(i) {
    cbind(i, adj_list[[i]])
  }))
  nug_from_edge_list(n, edges, ...)
}

#' Create an undirected graph from an adjacency matrix
#'
#' @param adj_mat A square symmetric matrix. Nonzero entries indicate edges;
#'   their values are used as edge weights. Self-loops are not allowed.
#' @param ... Additional arguments: `seed` (integer) for reproducible RNG.
#' @return A `NaturalUndirectedGraph` object.
#' @family graph-construction
#' @examples
#' # 4-vertex path graph from adjacency matrix
#' A <- matrix(0, 4, 4)
#' A[1, 2] <- A[2, 1] <- 1
#' A[2, 3] <- A[3, 2] <- 1
#' A[3, 4] <- A[4, 3] <- 1
#' g <- nug_from_adj_mat(A)
#' g$nedges()
#' @export
nug_from_adj_mat <- function(adj_mat, ...) {
  stopifnot(is.matrix(adj_mat))
  stopifnot(isSymmetric(adj_mat))
  if (any(diag(adj_mat) != 0)) {
    stop("Self-loops are not allowed")
  }
  n <- nrow(adj_mat)
  edges <- which(adj_mat != 0, arr.ind = TRUE)
  edges <- cbind(edges, adj_mat[edges])
  nug_from_edge_list(n, edges, ...)
}

#' Create a regular 2D grid graph
#'
#' Generates a grid graph with the specified dimensions and neighborhood order.
#' Order 1 gives rook adjacency (4-connected), order 2 gives queen adjacency
#' (8-connected).
#'
#' @param nrows Number of rows in the grid.
#' @param ncols Number of columns in the grid.
#' @param order Neighborhood order: `1` for rook (default) or `2` for queen.
#' @param ... Additional arguments: `seed` (integer) for reproducible RNG.
#' @return A `NaturalUndirectedGraph` object.
#' @family graph-construction
#' @examples
#' # 4x4 rook grid
#' g <- nug_from_grid(4, 4)
#' g$nvertices()
#' g$nedges()
#' @export
nug_from_grid <- function(nrows, ncols, order = 1L, ...) {
  nrows <- as.integer(nrows)
  ncols <- as.integer(ncols)
  order <- as.integer(order)
  graph_ptr <- nug_generate_regular_cpp(nrows, ncols, order)
  NaturalUndirectedGraph$new(.graph_ptr = graph_ptr, ...)
}
