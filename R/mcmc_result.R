#' MCMC Result
#'
#' R6 class containing posterior samples from [mcmc()]. Access samples via
#' the accessor methods below.
#'
#' @importFrom R6 R6Class
#' @export
MdgmResult <- R6::R6Class(
  classname = "MdgmResult",
  cloneable = FALSE,
  public = list(
    #' @description Create a new MdgmResult. Use [mcmc()] instead.
    #' @param raw List of raw MCMC output from C++.
    initialize = function(raw) {
      private$.raw <- raw
    },

    #' @description Get the latent field samples.
    #' @param iteration Optional integer iteration to extract (1-indexed).
    #'   If `NULL`, returns the full `n x J` matrix.
    #' @return Integer matrix (`n x J`) or integer vector (length `n`).
    z = function(iteration = NULL) {
      if (is.null(iteration)) return(private$.raw$z)
      private$.raw$z[, iteration]
    },

    #' @description Get the psi (dependence parameter) samples.
    #' @return Numeric vector of length `J`.
    psi = function() {
      private$.raw$psi
    },

    #' @description Get the emission parameter samples.
    #' @param iteration Optional integer iteration to extract (1-indexed).
    #'   If `NULL`, returns the full `n_colors x J` matrix.
    #' @return Numeric matrix (`n_colors x J`), numeric vector, or `NULL`
    #'   if standalone model.
    eta = function(iteration = NULL) {
      if (is.null(private$.raw$eta)) return(NULL)
      if (is.null(iteration)) return(private$.raw$eta)
      private$.raw$eta[, iteration]
    },

    #' @description Get the DAG parent vector samples.
    #' @param iteration Optional integer iteration to extract (1-indexed).
    #'   If `NULL`, returns the full `n x J` matrix.
    #' @return Integer matrix (`n x J`) or integer vector (length `n`).
    #'   Values are 1-indexed parent vertex IDs; `NA` indicates a root.
    dag = function(iteration = NULL) {
      if (is.null(iteration)) return(private$.raw$dag)
      private$.raw$dag[, iteration]
    },

    #' @description Compute edge inclusion probabilities.
    #' @param nug A `NaturalUndirectedGraph` object (the same one used
    #'   to create the model).
    #' @param burnin Number of initial iterations to discard (default 0).
    #' @return A data frame with columns `vertex1`, `vertex2`, and `prob`.
    edge_inclusion_probs = function(nug, burnin = 0L) {
      dag_mat <- private$.raw$dag
      n <- nrow(dag_mat)
      J <- ncol(dag_mat)
      start <- as.integer(burnin) + 1L
      if (start > J) stop("burnin exceeds number of iterations")
      dag_sub <- dag_mat[, start:J, drop = FALSE]
      J_eff <- ncol(dag_sub)

      edges <- list()
      for (v in seq_len(n)) {
        nbrs <- nug$neighbors(v)
        for (u in nbrs) {
          if (u > v) {
            edges[[length(edges) + 1L]] <- c(v, u)
          }
        }
      }

      v1 <- integer(length(edges))
      v2 <- integer(length(edges))
      prob <- numeric(length(edges))
      for (idx in seq_along(edges)) {
        e <- edges[[idx]]
        count <- sum(
          dag_sub[e[1], ] == e[2] | dag_sub[e[2], ] == e[1],
          na.rm = TRUE
        )
        v1[idx] <- e[1]
        v2[idx] <- e[2]
        prob[idx] <- count / J_eff
      }

      data.frame(vertex1 = v1, vertex2 = v2, prob = prob)
    },

    #' @description Get acceptance rates for MH steps.
    #' @return A named numeric vector with `psi` and `graph` rates.
    acceptance_rates = function() {
      J <- private$.raw$n_iterations - 1L
      c(
        psi = private$.raw$psi_accepted / J,
        graph = private$.raw$graph_accepted / J
      )
    },

    #' @description Print a summary of the MCMC results.
    #' @return Invisible `self`.
    summary = function() {
      J <- private$.raw$n_iterations
      n <- private$.raw$n_vertices
      nc <- private$.raw$n_colors
      rates <- self$acceptance_rates()
      cat(sprintf("MDGM MCMC Results\n"))
      cat(sprintf("  Vertices: %d, Colors: %d\n", n, nc))
      cat(sprintf("  Iterations: %d\n", J))
      cat(sprintf("  Psi acceptance rate: %.3f\n", rates["psi"]))
      cat(sprintf(
        "  Psi posterior mean: %.4f\n",
        mean(private$.raw$psi)
      ))
      if (!is.null(private$.raw$eta)) {
        eta_means <- rowMeans(private$.raw$eta)
        cat(sprintf(
          "  Eta posterior means: %s\n",
          paste(round(eta_means, 4), collapse = ", ")
        ))
      }
      invisible(self)
    }
  ),
  private = list(
    .raw = NULL
  )
)
