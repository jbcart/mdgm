# Suppress R CMD check notes for ggplot2 aes() variables
utils::globalVariables(c("iteration", "value", "parameter", "x", "y", "z"))

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
    #' @param emission_type Character string or NULL.
    #' @param nug A `NaturalUndirectedGraph` or NULL.
    initialize = function(raw, emission_type = NULL, nug = NULL) {
      private$.raw <- raw
      private$.emission_type <- emission_type
      private$.nug <- nug
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

    #' @description Get the emission parameters.
    #' @param iteration Optional integer iteration to extract (1-indexed).
    #' @return Named list appropriate to emission type, or `NULL` for
    #'   standalone models.
    emission_params = function(iteration = NULL) {
      if (is.null(private$.raw$theta)) return(NULL)
      et <- private$.emission_type
      if (is.null(et)) et <- "bernoulli"  # fallback

      theta_mat <- private$.raw$theta
      if (!is.null(iteration)) {
        theta_mat <- theta_mat[, iteration, drop = FALSE]
      }
      nc <- nrow(theta_mat)

      switch(et,
        bernoulli = {
          rownames(theta_mat) <- paste0("p_", seq_len(nc))
          list(p = theta_mat)
        },
        gaussian = {
          # theta stores mu (first half) and sigma^2 (second half)
          half <- nc %/% 2L
          mu_mat <- theta_mat[seq_len(half), , drop = FALSE]
          sigma2_mat <- theta_mat[half + seq_len(half), , drop = FALSE]
          rownames(mu_mat) <- paste0("mu_", seq_len(half))
          rownames(sigma2_mat) <- paste0("sigma2_", seq_len(half))
          list(mu = mu_mat, sigma2 = sigma2_mat)
        },
        poisson = {
          rownames(theta_mat) <- paste0("lambda_", seq_len(nc))
          list(lambda = theta_mat)
        },
        {
          rownames(theta_mat) <- paste0("theta_", seq_len(nc))
          list(theta = theta_mat)
        }
      )
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
    #' @param nug A `NaturalUndirectedGraph` object. If `NULL`, uses the
    #'   graph stored at construction time.
    #' @param burnin Number of initial iterations to discard (default 0).
    #' @return A data frame with columns `vertex1`, `vertex2`, and `prob`.
    edge_inclusion_probs = function(nug = NULL, burnin = 0L) {
      if (is.null(nug)) nug <- private$.nug
      if (is.null(nug)) stop("No graph available. Pass nug argument.")

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

    #' @description Compute MCMC diagnostics.
    #' @param burnin Number of initial iterations to discard (default 0).
    #' @return A named list with `rhat` and `ess` (effective sample size)
    #'   for psi and emission parameters.
    diagnostics = function(burnin = 0L) {
      start <- as.integer(burnin) + 1L
      J <- private$.raw$n_iterations
      if (start >= J) stop("burnin exceeds number of iterations")

      psi_chain <- private$.raw$psi[start:J]

      result <- list(
        psi = list(
          rhat = split_rhat(psi_chain),
          ess = ess(psi_chain)
        )
      )

      if (!is.null(private$.raw$theta)) {
        nc <- nrow(private$.raw$theta)
        theta_rhat <- numeric(nc)
        theta_ess <- numeric(nc)
        for (k in seq_len(nc)) {
          chain_k <- private$.raw$theta[k, start:J]
          theta_rhat[k] <- split_rhat(chain_k)
          theta_ess[k] <- ess(chain_k)
        }
        et <- if (!is.null(private$.emission_type)) private$.emission_type else "theta"
        param_name <- switch(et,
          bernoulli = "p",
          poisson = "lambda",
          gaussian = "mu_sigma2",
          "theta"
        )
        result[[param_name]] <- list(rhat = theta_rhat, ess = theta_ess)
      }

      result
    },

    #' @description Print a summary of the MCMC results.
    #' @param burnin Number of initial iterations to discard (default 0).
    #' @return Invisible named list with summary statistics.
    summary = function(burnin = 0L) {
      J <- private$.raw$n_iterations
      n <- private$.raw$n_vertices
      nc <- private$.raw$n_colors
      start <- as.integer(burnin) + 1L
      rates <- self$acceptance_rates()

      psi_post <- private$.raw$psi[start:J]

      cat(sprintf("MDGM MCMC Results\n"))
      cat(sprintf("  Vertices: %d, Colors: %d\n", n, nc))
      cat(sprintf("  Iterations: %d (burnin: %d)\n", J, burnin))
      cat(sprintf("  Psi acceptance rate: %.3f\n", rates["psi"]))
      cat(sprintf("  Psi posterior mean: %.4f (sd: %.4f)\n",
                  mean(psi_post), sd(psi_post)))

      out <- list(
        n_iter = J,
        n_vertices = n,
        n_colors = nc,
        burnin = burnin,
        acceptance_rates = rates,
        psi_mean = mean(psi_post),
        psi_sd = sd(psi_post)
      )

      if (!is.null(private$.raw$theta)) {
        theta_post <- private$.raw$theta[, start:J, drop = FALSE]
        theta_means <- rowMeans(theta_post)
        theta_sds <- apply(theta_post, 1, sd)
        et <- if (!is.null(private$.emission_type)) private$.emission_type else "unknown"
        n_theta <- nrow(theta_post)
        param_names <- switch(et,
          bernoulli = paste0("p_", seq_len(n_theta)),
          poisson = paste0("lambda_", seq_len(n_theta)),
          gaussian = c(paste0("mu_", seq_len(nc)),
                       paste0("sigma2_", seq_len(nc))),
          paste0("theta_", seq_len(n_theta))
        )
        cat(sprintf("  Emission type: %s\n", et))
        for (i in seq_along(theta_means)) {
          cat(sprintf("  %s posterior mean: %.4f (sd: %.4f)\n",
                      param_names[i], theta_means[i], theta_sds[i]))
        }
        out$emission_type <- et
        out$emission_param_means <- setNames(theta_means, param_names)
        out$emission_param_sds <- setNames(theta_sds, param_names)
      }

      # Diagnostics
      diag <- self$diagnostics(burnin = burnin)
      cat(sprintf("  Diagnostics:\n"))
      cat(sprintf("    psi — R-hat: %.4f, ESS: %.0f\n",
                  diag$psi$rhat, diag$psi$ess))

      if (!is.null(private$.raw$theta)) {
        # Find the theta diagnostics (keyed by emission type)
        diag_keys <- setdiff(names(diag), "psi")
        if (length(diag_keys) > 0) {
          theta_diag <- diag[[diag_keys[1]]]
          for (i in seq_along(param_names)) {
            cat(sprintf("    %s — R-hat: %.4f, ESS: %.0f\n",
                        param_names[i], theta_diag$rhat[i], theta_diag$ess[i]))
          }
        }
      }

      out$diagnostics <- diag

      invisible(out)
    },

    #' @description Plot MCMC trace plots and diagnostics.
    #' @param burnin Number of initial iterations to discard (default 0).
    #' @param which Character vector of what to plot: `"trace"` for trace
    #'   plots of psi and emission params, `"edge_inclusion"` for edge
    #'   inclusion probabilities (requires stored graph and igraph),
    #'   `"posterior_field"` for posterior mode of the latent field
    #'   (requires stored graph).
    #' @return Invisible list of ggplot/igraph plot objects.
    plot = function(burnin = 0L, which = "trace") {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 is required for plot(). Install it with install.packages('ggplot2').")
      }
      which <- match.arg(which, c("trace", "edge_inclusion", "posterior_field"),
                         several.ok = TRUE)
      start <- as.integer(burnin) + 1L
      J <- private$.raw$n_iterations

      plots <- list()

      if ("trace" %in% which) {
        # Psi trace plot
        psi_df <- data.frame(
          iteration = start:J,
          value = private$.raw$psi[start:J]
        )
        p_psi <- ggplot2::ggplot(psi_df, ggplot2::aes(x = iteration, y = value)) +
          ggplot2::geom_line(alpha = 0.6) +
          ggplot2::labs(title = "Trace: psi", x = "Iteration", y = "psi") +
          ggplot2::theme_minimal()
        print(p_psi)
        plots$psi <- p_psi

        # Emission params trace plots
        if (!is.null(private$.raw$theta)) {
          nc <- nrow(private$.raw$theta)
          et <- if (!is.null(private$.emission_type)) private$.emission_type else "theta"
          param_names <- switch(et,
            bernoulli = paste0("p_", seq_len(nc)),
            poisson = paste0("lambda_", seq_len(nc)),
            gaussian = c(paste0("mu_", seq_len(nc %/% 2L)),
                         paste0("sigma2_", seq_len(nc %/% 2L))),
            paste0("theta_", seq_len(nc))
          )
          theta_post <- private$.raw$theta[, start:J, drop = FALSE]
          iters <- start:J
          rows <- lapply(seq_len(nc), function(k) {
            data.frame(
              iteration = iters,
              value = theta_post[k, ],
              parameter = param_names[k]
            )
          })
          theta_df <- do.call(rbind, rows)
          p_theta <- ggplot2::ggplot(theta_df,
              ggplot2::aes(x = iteration, y = value, color = parameter)) +
            ggplot2::geom_line(alpha = 0.6) +
            ggplot2::labs(title = paste("Trace:", et, "parameters"),
                         x = "Iteration", y = "Value") +
            ggplot2::theme_minimal()
          print(p_theta)
          plots$emission <- p_theta
        }
      }

      if ("edge_inclusion" %in% which) {
        nug <- private$.nug
        if (is.null(nug)) {
          message("No graph stored; cannot plot edge inclusion probabilities.")
        } else if (!requireNamespace("igraph", quietly = TRUE)) {
          message("igraph is required for edge inclusion plot.")
        } else {
          eip <- self$edge_inclusion_probs(burnin = burnin)
          n_v <- nrow(private$.raw$z)

          # Build igraph object
          el <- as.matrix(eip[, c("vertex1", "vertex2")])
          g <- igraph::graph_from_edgelist(el, directed = FALSE)

          # Grid layout if square
          nside <- as.integer(sqrt(n_v))
          if (nside * nside == n_v) {
            coords <- cbind((seq_len(n_v) - 1) %% nside + 1,
                            nside - (seq_len(n_v) - 1) %/% nside)
          } else {
            coords <- igraph::layout_nicely(g)
          }

          # Color vertices by posterior mode of z
          z_post <- private$.raw$z[, (as.integer(burnin) + 1L):ncol(private$.raw$z),
                                    drop = FALSE]
          nc <- private$.raw$n_colors
          z_mode <- apply(z_post, 1, function(row) {
            tbl <- tabulate(row + 1L, nbins = nc)
            which.max(tbl) - 1L
          })
          pal <- if (nc == 2) c("#440154", "#fde725") else
                   grDevices::hcl.colors(nc, palette = "viridis")
          vcol <- pal[z_mode + 1L]

          igraph::plot.igraph(g, layout = coords,
               vertex.size = max(5, 30 / sqrt(n_v) * 2),
               vertex.label = NA,
               vertex.color = vcol,
               edge.width = eip$prob * 8,
               edge.color = "grey30",
               main = "Edge inclusion probabilities")
          plots$edge_inclusion <- eip
        }
      }

      if ("posterior_field" %in% which) {
        nug <- private$.nug
        if (is.null(nug)) {
          message("No graph stored; cannot plot posterior field.")
        } else {
          z_post <- private$.raw$z[, start:J, drop = FALSE]
          nc <- private$.raw$n_colors
          # Posterior mode per vertex
          z_mode <- apply(z_post, 1, function(row) {
            tbl <- tabulate(row + 1L, nbins = nc)
            which.max(tbl) - 1L
          })
          n <- length(z_mode)
          nside <- as.integer(sqrt(n))
          if (nside * nside == n) {
            field_df <- data.frame(
              x = rep(seq_len(nside), nside),
              y = rep(seq_len(nside), each = nside),
              z = factor(z_mode)
            )
            p_field <- ggplot2::ggplot(field_df,
                ggplot2::aes(x = x, y = y, fill = z)) +
              ggplot2::geom_tile() +
              ggplot2::scale_fill_viridis_d() +
              ggplot2::coord_equal() +
              ggplot2::labs(title = "Posterior mode of latent field",
                           fill = "Color") +
              ggplot2::theme_minimal()
            print(p_field)
            plots$field <- p_field
          } else {
            message("Posterior field plot only supported for square grids.")
          }
        }
      }

      invisible(plots)
    }
  ),
  private = list(
    .raw = NULL,
    .emission_type = NULL,
    .nug = NULL
  )
)

# --- Internal MCMC diagnostic helpers ---

#' Split R-hat for a single chain
#' @param x Numeric vector (single MCMC chain).
#' @return Scalar R-hat value.
#' @keywords internal
split_rhat <- function(x) {
  n <- length(x)
  if (n < 4L) return(NA_real_)
  half <- n %/% 2L
  chains <- list(x[seq_len(half)], x[(n - half + 1L):n])
  m <- 2L
  chain_n <- half

  chain_means <- vapply(chains, mean, double(1))
  chain_vars <- vapply(chains, stats::var, double(1))
  overall_mean <- mean(chain_means)

  B <- chain_n / (m - 1L) * sum((chain_means - overall_mean)^2)
  W <- mean(chain_vars)

  if (is.na(W) || W == 0) return(NA_real_)
  var_hat <- (chain_n - 1L) / chain_n * W + B / chain_n
  sqrt(var_hat / W)
}

#' Effective sample size via autocorrelation
#' @param x Numeric vector (single MCMC chain).
#' @return Scalar ESS estimate.
#' @keywords internal
ess <- function(x) {
  n <- length(x)
  if (n < 4L) return(NA_real_)
  ac <- stats::acf(x, lag.max = n - 1L, plot = FALSE)$acf[, , 1]
  # Geyer's initial positive sequence: sum consecutive pairs
  tau <- 1.0
  max_lag <- length(ac) - 1L
  t <- 1L
  while (t + 1L <= max_lag) {
    rho_pair <- ac[t + 1L] + ac[t + 2L]
    if (rho_pair < 0) break
    tau <- tau + 2.0 * rho_pair
    t <- t + 2L
  }
  max(1, n / tau)
}
