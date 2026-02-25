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
      if (is.null(private$.raw$eta)) return(NULL)
      et <- private$.emission_type
      if (is.null(et)) et <- "bernoulli"  # fallback

      eta_mat <- private$.raw$eta
      if (!is.null(iteration)) {
        eta_mat <- eta_mat[, iteration, drop = FALSE]
      }
      nc <- nrow(eta_mat)
      colnames_eta <- paste0("eta_", seq_len(nc))
      rownames(eta_mat) <- colnames_eta

      switch(et,
        bernoulli = list(eta = eta_mat),
        gaussian = {
          # eta stores mu (first nc rows) and sigma (next nc rows)
          half <- nc %/% 2L
          list(
            mu = eta_mat[seq_len(half), , drop = FALSE],
            sigma = eta_mat[half + seq_len(half), , drop = FALSE]
          )
        },
        poisson = {
          rownames(eta_mat) <- paste0("lambda_", seq_len(nc))
          list(lambda = eta_mat)
        },
        list(eta = eta_mat)
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

      if (!is.null(private$.raw$eta)) {
        nc <- nrow(private$.raw$eta)
        eta_rhat <- numeric(nc)
        eta_ess <- numeric(nc)
        for (k in seq_len(nc)) {
          chain_k <- private$.raw$eta[k, start:J]
          eta_rhat[k] <- split_rhat(chain_k)
          eta_ess[k] <- ess(chain_k)
        }
        et <- if (!is.null(private$.emission_type)) private$.emission_type else "eta"
        param_name <- switch(et,
          bernoulli = "eta",
          poisson = "lambda",
          gaussian = "mu_sigma",
          "eta"
        )
        result[[param_name]] <- list(rhat = eta_rhat, ess = eta_ess)
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

      if (!is.null(private$.raw$eta)) {
        eta_post <- private$.raw$eta[, start:J, drop = FALSE]
        eta_means <- rowMeans(eta_post)
        eta_sds <- apply(eta_post, 1, sd)
        et <- if (!is.null(private$.emission_type)) private$.emission_type else "unknown"
        param_names <- switch(et,
          bernoulli = paste0("eta_", seq_len(nc)),
          poisson = paste0("lambda_", seq_len(nc)),
          gaussian = c(paste0("mu_", seq_len(nc %/% 2L)),
                       paste0("sigma_", seq_len(nc %/% 2L))),
          paste0("eta_", seq_len(nc))
        )
        cat(sprintf("  Emission type: %s\n", et))
        for (i in seq_along(eta_means)) {
          cat(sprintf("  %s posterior mean: %.4f (sd: %.4f)\n",
                      param_names[i], eta_means[i], eta_sds[i]))
        }
        out$emission_type <- et
        out$emission_param_means <- setNames(eta_means, param_names)
        out$emission_param_sds <- setNames(eta_sds, param_names)
      }

      # Diagnostics
      diag <- self$diagnostics(burnin = burnin)
      cat(sprintf("  Psi R-hat: %.4f, ESS: %.0f\n",
                  diag$psi$rhat, diag$psi$ess))
      out$diagnostics <- diag

      invisible(out)
    },

    #' @description Plot MCMC trace plots and diagnostics.
    #' @param burnin Number of initial iterations to discard (default 0).
    #' @param which Character vector of what to plot: `"trace"` for trace
    #'   plots of psi and emission params, `"posterior_field"` for posterior
    #'   mean of the latent field (requires stored graph).
    #' @return Invisible `NULL`.
    plot = function(burnin = 0L, which = "trace") {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 is required for plot(). Install it with install.packages('ggplot2').")
      }
      which <- match.arg(which, c("trace", "posterior_field"), several.ok = TRUE)
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
        if (!is.null(private$.raw$eta)) {
          nc <- nrow(private$.raw$eta)
          et <- if (!is.null(private$.emission_type)) private$.emission_type else "eta"
          param_names <- switch(et,
            bernoulli = paste0("eta_", seq_len(nc)),
            poisson = paste0("lambda_", seq_len(nc)),
            paste0("eta_", seq_len(nc))
          )
          eta_post <- private$.raw$eta[, start:J, drop = FALSE]
          iters <- start:J
          rows <- lapply(seq_len(nc), function(k) {
            data.frame(
              iteration = iters,
              value = eta_post[k, ],
              parameter = param_names[k]
            )
          })
          eta_df <- do.call(rbind, rows)
          p_eta <- ggplot2::ggplot(eta_df,
              ggplot2::aes(x = iteration, y = value, color = parameter)) +
            ggplot2::geom_line(alpha = 0.6) +
            ggplot2::labs(title = paste("Trace:", et, "parameters"),
                         x = "Iteration", y = "Value") +
            ggplot2::theme_minimal()
          print(p_eta)
          plots$emission <- p_eta
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

  if (W == 0) return(NA_real_)
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
