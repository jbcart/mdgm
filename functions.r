# libraries
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(sf))
sf_use_s2(FALSE)
suppressPackageStartupMessages(library(progress))
suppressPackageStartupMessages(library(truncdist))
suppressPackageStartupMessages(library(GiRaF))
suppressPackageStartupMessages(library(igraph))

# functions for finding neighbors and matrices
st_rook <- function(a, b=a) st_relate(a, b, pattern="F***1****")
st_queen <- function(a, b=a) st_relate(a, b, pattern="F***T****")
nbh_rook <- function(a) 1*as.matrix(st_rook(a))
nbh_queen <- function(a) 1*as.matrix(st_queen(a))

graph_from_nug <- function(nugmat) {
  nugmat[upper.tri(nugmat)] <- 0
  nbh_partial <- apply(nugmat, 1, function(x) (1:n)[x==1])
  el <- Reduce(rbind, lapply(1:length(nbh_partial), function(i) {
    if(length(nbh_partial[[i]]) == 0) matrix(,nrow=0,ncol=2)
    else cbind(nbh_partial[[i]],i)
  }))
  graph_from_edgelist(el, directed=FALSE)
}

graph_from_dag <- function(dag) {
  el <- Reduce(rbind, lapply(1:length(dag), function(i) {
    if(length(dag[[i]]) == 0) matrix(,nrow=0,ncol=2)
    else cbind(dag[[i]],i)
  }))
  graph_from_edgelist(el, directed=TRUE)
}

# function for weights for rootd dag algorithm
make_weights_rdag <- function(nug_rook, nug_queen) {
  if (missing(nug_queen)) {
    lapply(nug_rook, function(x) 1+numeric(length(x)))
  } else {
    weights <- lapply(nug_queen, function(x) 1+numeric(length(x)))
    for (i in 1:length(weights)) {
      weights[[i]][!(nug_queen[[i]] %in% nug_rook[[i]])] <- 2
    }
    weights
  }
}

# generate list representation of dependence structure
# this assumes lattice is labeled from bottom left corner, 
# labels increases left to right.
nbh_mrf_q <- function(nrow, ncol) {
  get_index <- function(i) {
    l <- i %% ncol != 1
    r <- i %% ncol != 0
    t <- (i-1) %/% ncol != (nrow-1)
    b <- (i-1) %/% ncol != 0
    i_t <- (i %/% ncol + 1) * ncol + i %% ncol
    i_b <- (i %/% ncol - 1) * ncol + i %% ncol
    neighbors <- c((i_t-1)*l*t, i_t*t, (i_t+1)*r*t,
                   (i-1)*l, (i+1)*r,
                   (i_b-1)*l*b, i_b*b, (i_b+1)*r*b)
    neighbors[neighbors != 0]
  }
  lapply(1:(ncol*nrow), get_index)
}

nbh_mrf_r <- function(nrow, ncol) {
  get_index <- function(i) {
    l <- i %% ncol != 1
    r <- i %% ncol != 0
    t <- (i-1) %/% ncol != (nrow-1)
    b <- (i-1) %/% ncol != 0
    i_t <- (i %/% ncol + 1) * ncol + i %% ncol
    i_b <- (i %/% ncol - 1) * ncol + i %% ncol
    neighbors <- c(i_t*t,
                   (i-1)*l, (i+1)*r,
                   i_b*b
                  )
    neighbors[neighbors != 0]
  }
  lapply(1:(ncol*nrow), get_index)
}

# Generates a rooted DAG that is "faithful" 
#  to the conditional independence assumptions
root_rq <- function(root, nrow, ncol) {
  if (root %% nrow == 0) c(root %/% nrow, nrow)
  else c(root %/% nrow + 1, root %% nrow)
}

generate_rdag <- function(nug, weights, root) {
  n <- length(nug)
  # create labels from root
  labels <- numeric(n)
  labeled <- logical(n)
  labeled[root] <- TRUE
  next_vertices <- root
  while (!all(labeled)) {
    connections <- nug[next_vertices]
    df <- data.frame(v=Reduce(c, connections),
                     d=Reduce(c, weights[next_vertices]) + 
                       rep(labels[next_vertices], times=sapply(connections, length)))
    df <- df[!labeled[df$v],]
    df <- df[order(df$v,df$d),]
    df <- df[!duplicated(df$v),]
    labels[df$v] <- df$d
    labeled[df$v] <- TRUE
    next_vertices <- df$v
  }
  list(dag=lapply(1:length(nug), function(i) nug[[i]][labels[nug[[i]]] < labels[i]]), labels=labels)
}

nbh_mdgm_r <- function(root, nrow, ncol) {
  n <- ncol*nrow
  root_id <- root_rq(root, nrow, ncol)
  r_root <- root_id[1]
  q_root <- root_id[2]
  get_index_gc <- function(i) {
    id <- root_rq(i, nrow, ncol)
    r <- id[1]
    q <- id[2]
    i_t <- (i %/% ncol + 1) * ncol + i %% ncol
    i_b <- (i %/% ncol - 1) * ncol + i %% ncol
    if (q < q_root) {
      if (r > r_root) c(i+1, i_b, i_b+1)
      else if (r == r_root) i+1
      else c(i_t, i_t+1, i+1)
    } else if (q == q_root) {
      if (r > r_root) i_b
      else if (r == r_root) numeric(0)
      else i_t
    } else if (q > q_root) {
      if (r > r_root) c(i-1, i_b-1, i_b)
      else if (r == r_root) i-1
      else c(i-1, i_t-1, i_t)
    }
  }
  lapply(1:n, get_index_gc)
}

# generate random uniform spanning tree
nbh_rust <- function(nbh) {
  n <- length(nbh)
  in_tree <- logical(n)
  root <- sample(1:n, 1)
  in_tree[root] <- TRUE
  spanning_tree <- lapply(1:n, function(i) numeric())

  while (sum(in_tree) < n) {
    starters <- which(!in_tree)
    random_walk <- ifelse(length(starters) == 1, starters, sample(starters, 1))
    next_node <- sample(nbh[[random_walk]], 1)
    while (!in_tree[next_node]) {
      if (next_node %in% random_walk) {
        loop_start <- match(next_node, random_walk)
        random_walk <- random_walk[1:loop_start]
      } else {
        random_walk <- c(random_walk, next_node)
      }
      next_node <- sample(nbh[[next_node]], 1)
    }
    random_walk <- c(random_walk, next_node)
    for (i in 2:length(random_walk)) {
      spanning_tree[[random_walk[i-1]]] <- random_walk[i]
    }
    in_tree[random_walk] <- TRUE
  }
  spanning_tree
}

# sample a spanning tree from posterior distribution given z and psi
nbh_st_post <- function(nbh, z, psi) {
  n <- length(nbh)
  in_tree <- logical(n)
  root <- sample(1:n, 1)
  in_tree[root] <- TRUE
  spanning_tree <- lapply(1:n, function(i) numeric())

  while (sum(in_tree) < n) {
    starters <- which(!in_tree)
    random_walk <- ifelse(length(starters) == 1, starters, sample(starters, 1))
    prob_next_node <- exp(psi * (z[random_walk] == z[nbh[[random_walk]]]))
    if (length(nbh[[random_walk]]) == 1){
      next_node <- nbh[[random_walk]]
    } else {
      next_node <- sample(nbh[[random_walk]], 1, prob=prob_next_node)
    }
    while (!in_tree[next_node]) {
      if (next_node %in% random_walk) {
        loop_start <- match(next_node, random_walk)
        random_walk <- random_walk[1:loop_start]
      } else {
        random_walk <- c(random_walk, next_node)
      }
      prob_next_node <- exp(psi * (z[next_node] == z[nbh[[next_node]]]))
      if (length(nbh[[next_node]]) == 1){
        next_node <- nbh[[next_node]]
      } else {
        next_node <- sample(nbh[[next_node]], 1, prob=prob_next_node)
      }
    }
    random_walk <- c(random_walk, next_node)
    for (i in 2:length(random_walk)) {
      spanning_tree[[random_walk[i-1]]] <- random_walk[i]
    }
    in_tree[random_walk] <- TRUE
  }
  spanning_tree
}

# generate ao mdgm
nbh_mdgm_ao <- function(nbh, perm) {
  n <- length(nbh)
  if (missing(perm)) perm <- sample(1:n)
  nbh_out <- lapply(1:n, function(i) numeric())
  for (i in 1:n) {
    j <- perm[i]
    nbh_out[[j]] <- nbh[[j]][nbh[[j]] %in% perm[1:i]]
  }
  nbh_out
}

# simulate color conditional on parents/(nieghbors)
sim_cond_pa <- function(alpha, psi, pa, n_colors) {
  cat <- 0:(n_colors-1)
  if (length(pa) == 0) {
    w <- exp(alpha)
  } else {
    w <- exp(psi*table(factor(pa, levels=cat)) + alpha)
  }
  sample(cat, 1, prob=w)
}

# simulate color conditional on the parents given a uniform draw
# used for CFTP function
# currently implements for n_colors = 2 (Ising model)
sim_cond_pa_U <- function(alpha, psi, pa, u) {
  n_colors = 2
  cat <- 0:(n_colors-1)
  if (length(pa) == 0) {
    w <- exp(alpha)
  } else {
    w <- exp(psi*table(factor(pa, levels=cat)) + alpha)
  }
  ifelse(u < w[1]/sum(w), 0, 1)
}

# Coupling from the Past (CFTP) for exact sampling.
# currently implemented for the Ising model (k=2)
# nbh   - a neighborhood matrix
# alpha - marginal probability of each color
# psi   - inverse temperature (dependence parameter)
# init  - a list of two initialization points, if missing then 
#         colors are independently and randomly generated
# M     - Number of iterations to couple from the past
rmrf_cftp <- function(nbh, alpha = c(0, 0), psi=0.2, init, M=100, suppress_msg=FALSE) {
  n_colors <- length(alpha)
  if (is.matrix(nbh)) {
    n <- nrow(nbh)
    index <- 1:n
    nbh <- apply(nbh, 1, function(x) index[x == 1])
  } else {
    n <- length(nbh)
  }
  if (missing(init)) {
    init <- list()
    init[[1]] <- sapply(1:n, function(i) sim_cond_pa(alpha, psi, numeric(0), n_colors))
    init[[2]] <- sapply(1:n, function(i) sim_cond_pa(alpha, psi, numeric(0), n_colors))
  }
  color1 <- init[[1]]
  color2 <- init[[2]]
  U <- NULL
  B <- 0
  while (!all(color1 == color2)) {
    B <- B + M
    if (!suppress_msg) {
      cat(sprintf("MCMC chains have not coalesced!\nGoing back an additional %d steps for a total of %d steps.\n", M, B))
    }
    color1 <- init[[1]]
    color2 <- init[[2]]
    U <- cbind(matrix(runif(n*M), nrow=n, ncol=M), U)
    pb <- progress_bar$new(total=B)
    for (j in 1:B) {
      for (i in 1:n) {
        color1[i] <- sim_cond_pa_U(alpha, psi, color1[nbh[[i]]], U[i,j])
        color2[i] <- sim_cond_pa_U(alpha, psi, color2[nbh[[i]]], U[i,j])
      }
      pb$tick()
    }
  }
  list(color=color1)
}

# generate observed data
sample_emission <- function(z, eta, size=1) {
  p <- eta[z+1]
	fun <- function(m,s) {
		if (s!=0) {
		 	sample(c(0,1), size=s, replace=TRUE, prob=c(1-m, m))
		} else {
			NA
		}
	}
	mapply(fun, p, size, SIMPLIFY=FALSE)  
}

# Functions for MCMC
# probability of observed color conditional on parents/(neighbors)
d_cond_pa <- function(z_i, alpha, psi, pa) {
  n_colors <- length(alpha)
  cat <- 0:(n_colors-1)
  if (length(pa) == 0) {
    w <- exp(alpha)
  } else {
    w <- exp(psi*table(factor(pa, levels=cat)) + alpha)
  }
  w[z_i + 1]/sum(w)
}

# same as above but outputs probabilities for z[i] = {0,1}
p_cond_pa <- function(alpha, psi, pa) {
  n_colors <- length(alpha)
  cat <- 0:(n_colors-1)
  if (length(pa) == 0) {
    w <- exp(alpha)
  } else {
    w <- exp(psi*table(factor(pa, levels=cat)) + alpha)
  }
  w/sum(w)
}

# full conditionals for DAG are different than full conditionals for a MRF
p_full_cond_dag <- function(z, i, nbh, alpha, psi) {
  n_colors <- length(alpha)
  cat <- 0:(n_colors-1)
  gamma <- alpha[2]

  # denominator
  de <- which(sapply(nbh, function(x) i %in% x))
  if (length(de) == 0) {
    w_d <- c(1,1)
  } else {
    pa_of_de <- lapply(nbh[de], function(x) x[x!=i])
    w_d <- sapply(c(0,1), function(k) prod(
      sapply(pa_of_de, function(j) 
        exp(gamma + psi*(sum(z[j] == 1) + (k == 1))) +
        exp(psi*(sum(z[j] == 0) + (k == 0)))
    )))
  }

  # numerator
  nbs <- c(z[nbh[[i]]], z[de]) 
  if (length(nbs) == 0) {
    w <- exp(alpha)/w_d 
  } else {
    w <- exp(psi*table(factor(nbs, levels=cat)) + alpha)/w_d
  }
  if (any(is.na(w)) || length(w) != n_colors) stop("error with prob weights")
  
  if (any(is.na(w/sum(w)))) {
    w + 1
  } else {
    w/sum(w)
  }
}

# Numerator Likelihood for MRF (i.e. unnormalized)
# used for MH update of psi in MCMC with Moller algorithm
nl_mrf <- function(z, nbh, alpha, psi) {
  s <- sum(sapply(1:length(z), function(i) sum(z[i] == z[nbh[[i]]])))/2
  exp(sum(alpha[z+1]) + psi*s)
}

# Numerator log-Likelihood for MRF 
nll_mrf <- function(z, nbh, alpha, psi) {
  s <- sum(sapply(1:length(z), function(i) sum(z[i] == z[nbh[[i]]])))/2
  sum(alpha[z+1]) + psi*s
}

# Psuedo-Likelihood for mrf
# (True likelihood if nbh is a DAG)
pl_mrf <- function(z, nbh, alpha, psi) {
  cond_prob <- sapply(1:length(z),
  function(i) d_cond_pa(z[i],alpha,psi,z[nbh[[i]]]))
  prod(cond_prob)
}

# Psuedo-log-Likelihood for mrf
# (True likelihood if nbh is a DAG)
pll_mrf <- function(z, nbh, alpha, psi) {
  cond_prob <- sapply(1:length(z), function(i) d_cond_pa(z[i],alpha,psi,z[nbh[[i]]]))
  sum(log(cond_prob))
}

# likelihood of emissions
l_emission <- function(y, z, eta) {
  p <- eta[z+1]
	out <- Reduce(c, lapply(1:length(p), function(i) c(1-p[i], p[i])[y[[i]]+1]))
  prod(out, na.rm=TRUE)
}

# log likelihood of emissions
ll_emission <- function(y, z, eta) {
  p <- eta[z+1]
	out <- Reduce(c, lapply(1:length(p), function(i) c(1-p[i], p[i])[y[[i]]+1]))
  sum(log(out), na.rm=TRUE)
}

# priors
prior_psi <- function(psi) {
  dtrunc(psi, "cauchy", a=0)
}

prior_eta <- function(eta,a,b) {
  (eta[1] < eta[2]) * prod(dbeta(eta,a,b))
}

# log priors
log_prior_psi <- function(psi) {
  log(dtrunc(psi, "cauchy", a=0))
}

log_prior_eta <- function(eta,a,b) {
  log(eta[1] < eta[2]) + sum(dbeta(eta,a,b,log=TRUE))
}

# MCMC functions
# y_i|z_i,eta ~ bern(eta_{z[i]})
# eta ~ beta(eta[1],a,b) * beta(eta[2],a,b)
# z|psi ~ mrf(psi)
# or
# z|psi,dag ~ dagr(psi, dag)
# psi ~ gamma(c,d)
# 
# update latent variables (same for pl and complete mrf)
update_z_mrf <- function(y,z,nbh,psi,eta,alpha) {
  n <- length(z)
  for (i in sample(1:n)) {
		if (any(is.na(y[[i]]))) {
			w <- p_cond_pa(alpha,psi,z[nbh[[i]]])
		} else {
			w <- sapply(eta, function(p) prod(y[[i]]*p + (1-y[[i]])*(1-p))) * p_cond_pa(alpha,psi,z[nbh[[i]]])
		}
    z[i] <- sample(c(0,1),1,prob=w)
  }
  z
}

update_z_dag <- function(y,z,nbh,psi,eta,alpha) {
  n <- length(z)
  for (i in sample(1:n)) {
		if (any(is.na(y[[i]]))) {
			w <- p_full_cond_dag(z,i,nbh,alpha,psi)
		} else {
			w <- sapply(eta, function(p) prod(y[[i]]*p + (1-y[[i]])*(1-p))) * p_full_cond_dag(z,i,nbh,alpha,psi)
		}
    if (any(is.na(w)) || length(w) != length(alpha)) stop("error with probs")
    z[i] <- sample(c(0,1),1,prob=w)
  }
  z
}

# Update psi (metropolis step) 
# **same function for dag but neighborhood structure is different
update_psi_pl <- function(psi, psi_tune, z, nbh, alpha) {
  proposal <- psi + rnorm(1,0,sd=psi_tune)
  log_proposal <- pll_mrf(z, nbh, alpha, proposal) + log_prior_psi(proposal)
  log_current <- pll_mrf(z, nbh, alpha, psi) + log_prior_psi(psi)
  ifelse(log(runif(1)) < log_proposal - log_current, proposal, psi)
}

update_psi_mrf <- function(psi, psi_tune, z, nrow, ncol, nbh, alpha, M) {
  proposal <- psi + rnorm(1,0,sd=psi_tune)
  aux <- as.vector(exact.mrf(h=nrow, w=ncol, param=proposal, nei=8))
  log_proposal <- nll_mrf(z, nbh, alpha, proposal) + nll_mrf(aux, nbh, alpha, psi) + log_prior_psi(proposal)
  log_current <- nll_mrf(z, nbh, alpha, psi) + nll_mrf(aux, nbh, alpha, proposal) + log_prior_psi(psi)
  ifelse(log(runif(1)) < log_proposal - log_current, proposal, psi)
}

# update eta
# same for all schemes
update_eta <- function(y,z,eta,a,b) {
	n_z0 <- sum(!is.na(Reduce(c, y[z==0])))
	n_z1 <- sum(!is.na(Reduce(c, y[z==1])))
	data_z0 <- sum(Reduce(c, y[z==0]), na.rm=TRUE)
	data_z1 <- sum(Reduce(c, y[z==1]), na.rm=TRUE)
  eta[1] <- rtrunc(1, "beta", b=eta[2], shape1=a + data_z0, shape2=b + n_z0 - data_z0)
  eta[2] <- rtrunc(1, "beta", a=eta[1], shape1=a + data_z1, shape2=b + n_z1 - data_z1)
  eta
}

## Full MCMC
# y can be one of the following:
# - a list, where each element is the ith locations responses
# - a matrix, where each row is the ith locations observed responses
# - a vector, a single response at each location
# - a data frame with a column $y and $id, where $id is the location id
# Missing values should be indicated with NA
discrete_hmrf_mcmc <- function(y, nug, inits, psi_tune, M=100, J=1000, a=1, b=1,
 graph_prior=c("mrf_pl", "mrf_exact", "mdgm_ao", "mdgm_r", "mdgm_st")[1],
 graph_mh=c("independent", "random_walk")[1], weights=NULL, update_psi=TRUE, update_eta=TRUE, nrow=16, ncol=16) {
  # initialization
	if (is.numeric(y)) {
		y <- as.list(y)	
	} else if (is.matrix(y)) {
		y <- apply(y,1,c,simplify=FALSE)
	} else if (is.data.frame(y)) {
		y <- split(y, f=~id)
		y <- lapply(y, function(x) x$y)
	}
  n <- length(y)
  nbh_ref <- nbh <- nug
  Z <- matrix(0,nrow=n,ncol=J)
  Z[,1] <- inits$z
  eta <- matrix(0,2,ncol=J)
  eta[,1] <- inits$eta

  # alpha assumed to be zero, not updated
  alpha <- c(0,0)
  psi <- numeric(J)
  
  if (graph_prior %in% c("mdgm_r") & is.null(weights)) {
    stop("list of edge weights must be supplied for graph_prior=mdgm_r")
  }

  if (graph_prior %in% c("mdgm_r")) {
    root <- numeric(J)
    root[1] <- sample(1:n, 1)
  } else root <- numeric()

  if (graph_prior %in% c("mdgm_ao")) {
    perm <- matrix(0,nrow=n,ncol=J)
  } else perm <- numeric()

  if (graph_prior %in% c("mdgm_st")) {
    st <- matrix(0,nrow=n,ncol=J)
  } else st <- numeric()

  psi[1] <- inits$psi
  pb <- progress_bar$new(total=J)
  acc_psi <- 0
  acc_graph <- 0

  for (j in 2:J) {
    # For PL and MRF - graph does not change, just update Z.
    # For mdgm priors, update the graph (random walk or independent MH step)
    # then update Z
    if (graph_prior %in% c("mrf_pl", "mrf_exact")) {
      Z[,j] <- update_z_mrf(y, Z[,j-1], nbh, psi[j-1], eta[,j-1], alpha)
    } else if (graph_prior %in% c("mdgm_r", "mdgm_ao", "mdgm_st")) {
      if (graph_prior == "mdgm_r") {
        if (graph_mh == "independent") {
          root_prop <- sample(1:n, 1)
          ratio <- 0
        } else if (graph_mh == "random_walk") {
          root_prop <- sample(nbh_ref[[root[j-1]]], 1)
          ratio <- log(1/length(nbh_ref[[root_prop]])) - log(1/length(nbh_ref[[root[j-1]]]))
        }
        if (!is.list(nug)) {
          nbh_prop <- nbh_mdgm_r(root_prop, nrow, ncol)
        } else {
          nbh_prop <- generate_rdag(nbh_ref, weights, root_prop)$dag
        }
        ratio <- ratio + pll_mrf(Z[,j-1], nbh_prop, alpha, psi[j-1]) - pll_mrf(Z[,j-1], nbh, alpha, psi[j-1])
        if (log(runif(1)) < ratio) {
          nbh <- nbh_prop
          root[j] <- root_prop
          acc_graph <- acc_graph + 1
        } else {
          root[j] <- root[j-1]
        }
      } else if (graph_prior == "mdgm_ao") {
        # currently only independent MH step implemented
        perm_prop <- sample(1:n)
        nbh_prop <- nbh_mdgm_ao(nbh_ref, perm_prop)
        ratio <- pll_mrf(Z[,j-1], nbh_prop, alpha, psi[j-1]) - pll_mrf(Z[,j-1], nbh, alpha, psi[j-1])
        if (log(runif(1)) < ratio) {
         nbh <- nbh_prop
         perm[,j] <- perm_prop
         acc_graph <- acc_graph + 1
        } else {
         perm[,j] <- perm[,j-1]
        } 
      } else if (graph_prior == "mdgm_st") {
        nbh <- nbh_st_post(nbh_ref, Z[,j-1], psi[j-1])
        st[,j] <- append(Reduce(c, nbh), 0, match(0, sapply(nbh, length))-1)
      }
      Z[,j] <- update_z_dag(y, Z[,j-1], nbh, psi[j-1], eta[,j-1], alpha)
    } else {
      stop("Invalid graph_prior")
    }
    # update psi (same for psuedo-likeihood and dag)
    if (update_psi) {
      if (graph_prior == "mrf_exact") {
        proposal <- psi[j-1] + rnorm(1,0,sd=psi_tune)
				# exact.mrf has problems for nei=8 :( 
				# follow up - problems when param is above "citical value"
        aux <- as.vector(exact.mrf(h=nrow, w=ncol, param=proposal, nei=4))
        # wierd byproduct of exact.mrf function
        # aux[aux != 0 | aux != 1] <- 1
				# aux <- rmrf_cftp(nbh, psi=proposal, M=M)$color
        log_proposal <- nll_mrf(Z[,j], nbh, alpha, proposal) + nll_mrf(aux, nbh, alpha, psi[j-1]) + log_prior_psi(proposal)
        log_current <- nll_mrf(Z[,j], nbh, alpha, psi[j-1]) + nll_mrf(aux, nbh, alpha, proposal) + log_prior_psi(psi[j-1])
        psi[j] <- ifelse(log(runif(1)) < log_proposal - log_current, proposal, psi[j-1])
        # update_psi_mrf(psi[j-1], psi_tune, Z[,j], nrow, ncol, nbh, alpha, M)
      } else {
        psi[j] <- update_psi_pl(psi[j-1], psi_tune, Z[,j], nbh, alpha)
      }
      if (psi[j] != psi[j-1]) acc_psi <- acc_psi + 1
    } else {
      psi[j] <- psi[j-1]
    }
    # update eta
		if (update_eta) {
			eta[,j] <- update_eta(y, Z[,j], eta[,j-1], a, b)
		} else {
			eta[,j] <- eta[,j-1]
		}
		pb$tick()
  }
  list(Z=Z, eta=eta, psi=psi, acc_psi=acc_psi/J, acc_graph=acc_graph/J, root=root, perm=perm, st=st)
}

## Functions for evaluating performance
bbCount <- function(nbh, z) {
  sapply(seq_along(nbh), function(i) (z[nbh[[i]]] == 1) & (z[i] == 1)) %>%
    Reduce(sum,.)/2
}

bwCount <- function(nbh, z) {
  sapply(seq_along(nbh), function(i) z[nbh[[i]]] != z[i]) %>%
    Reduce(sum,.)/2
}

wwCount <- function(nbh, z) {
  sapply(seq_along(nbh), function(i) (z[nbh[[i]]] == 0) & (z[i] == 0)) %>%
    Reduce(sum,.)/2
}

matchesCount <- function(nbh, z) {
  sapply(seq_along(nbh), function(i) z[nbh[[i]]] == z[i]) %>%
    Reduce(sum,.)/2
}

edgesCountMRF <- function(nbh) {
  sum(sapply(nbh, length))/2
}

edgesCountDAG <- function(nbh) {
  sum(sapply(nbh, length))
}

moransI_cressie <- function(nbh, z) {
  sum(z * sapply(nbh, function(x) mean(z[x]))) / sum(z)
}

moransI <- function(nbh, z) {
  bbCount(nbh, z) * length(z) / (edgesCountMRF(nbh) * sum(z))
}

