#pragma once

#include <cmath>
#include <cstddef>
#include <mdgm/model.hpp>
#include <mdgm/observations.hpp>
#include <mdgm/rng.hpp>
#include <vector>

namespace mdgm {

// Half-Cauchy log-prior: log(2 / (pi * (1 + psi^2))) for psi > 0
inline double LogPriorPsi(double psi) {
  if (psi <= 0.0) return -1e300;
  return std::log(2.0 / M_PI) - std::log(1.0 + psi * psi);
}

struct McmcConfig {
  std::size_t n_iterations;
  double psi_tune;
  std::vector<double> emission_prior_params;
  bool store_z = false;
};

struct McmcSamples {
  std::vector<int> z;                    // n x J column-major (empty if !store_z)
  std::vector<double> psi;               // length J
  std::vector<double> theta;             // n_theta x J
  std::vector<std::size_t> dag_data;     // n x J compact DAG storage
  std::size_t psi_accepted;
  std::size_t graph_accepted;
  std::size_t n_iterations;
  std::size_t n_vertices;
  std::size_t n_colors;
  std::size_t n_theta;  // number of emission params per iteration

  // Always computed (cheap):
  std::vector<std::size_t> alloc;        // n x k allocation counts (post-burnin)
  std::vector<double> sufficient_stat;   // length J: T(z) same-color edge count
  std::vector<int> z_map;                // length n: joint MAP configuration
  double log_posterior_map = -1e300;      // score of z_map
  std::size_t map_iteration = 0;         // iteration index of z_map
};

McmcSamples RunMcmc(
    Model& model,
    const Observations& y,
    const std::vector<int>& z_init,
    double psi_init,
    const std::vector<double>& theta_init,
    const McmcConfig& config,
    RNG& rng);

}  // namespace mdgm
