#include <cmath>
#include <cstddef>
#include <mdgm/emission.hpp>
#include <mdgm/mcmc.hpp>
#include <mdgm/model.hpp>
#include <mdgm/observations.hpp>
#include <mdgm/rng.hpp>
#include <span>
#include <vector>

namespace mdgm {

McmcSamples RunMcmc(
    Model& model,
    const Observations& y,
    const std::vector<int>& z_init,
    double psi_init,
    const std::vector<double>& theta_init,
    const McmcConfig& config,
    RNG& rng) {
  const std::size_t n = model.nvertices();
  const std::size_t nc = model.ncolors();
  const std::size_t J = config.n_iterations;

  // Initialize theta
  std::vector<double> theta(theta_init);
  if (theta.empty() && model.has_emission()) {
    theta.resize(nc, 0.5);  // default if not provided
  }

  // Number of emission parameters per iteration
  const std::size_t n_theta = theta.size();

  // Allocate output storage
  McmcSamples samples;
  samples.n_iterations = J;
  samples.n_vertices = n;
  samples.n_colors = nc;
  samples.n_theta = n_theta;
  if (config.store_z) {
    samples.z.resize(n * J);
  }
  samples.psi.resize(J);
  samples.theta.resize(model.has_emission() ? n_theta * J : 0);
  samples.dag_data.resize(n * J);
  samples.psi_accepted = 0;
  samples.graph_accepted = 0;

  // Always-computed summaries
  samples.alloc.resize(n * nc, 0);
  samples.sufficient_stat.resize(J);
  samples.z_map.resize(n);
  samples.log_posterior_map = -1e300;
  samples.map_iteration = 0;

  // Initialize iteration 0
  std::vector<int> z(z_init);
  double psi = psi_init;

  if (config.store_z) {
    for (std::size_t i = 0; i < n; ++i) samples.z[i] = z[i];
  }
  samples.psi[0] = psi;
  for (std::size_t k = 0; k < n_theta; ++k) samples.theta[k] = theta[k];
  model.StoreDagSample(samples.dag_data, 0);

  // Accumulate summaries for iteration 0
  for (std::size_t i = 0; i < n; ++i) {
    samples.alloc[i * nc + static_cast<std::size_t>(z[i])] += 1;
  }
  samples.sufficient_stat[0] = model.SufficientStatistic(z);

  // Compute initial MAP score
  {
    double score = model.SpatialLogLikelihood(z, psi) + LogPriorPsi(psi);
    if (model.has_emission()) {
      score += EmissionLogLikelihood(y, z, theta, model.emission_type());
    }
    samples.log_posterior_map = score;
    samples.z_map.assign(z.begin(), z.end());
    samples.map_iteration = 0;
  }

  // MCMC iterations
  for (std::size_t j = 1; j < J; ++j) {
    // Step 1: Update graph structure
    model.UpdateGraph(z, psi, rng);
    model.StoreDagSample(samples.dag_data, j);

    // Step 2: Update z (Gibbs scan, only if hierarchical)
    if (model.has_emission()) {
      auto perm = rng.permutation(n);
      for (std::size_t idx = 0; idx < n; ++idx) {
        std::size_t i = perm[idx];
        auto probs = model.ZFullConditional(z, i, psi, y, theta);
        z[i] = static_cast<int>(rng.discrete(probs));
      }
    }

    // Step 3: Update psi (delegates to spatial field — default is MH random walk)
    psi = model.UpdatePsi(z, psi, config.psi_tune, samples.psi_accepted, rng);

    // Step 4: Update emission params (only if hierarchical)
    if (model.has_emission()) {
      theta = model.UpdateEmissionParams(y, z, theta,
                                          config.emission_prior_params, rng);
    }

    // Store full z (optional)
    if (config.store_z) {
      for (std::size_t i = 0; i < n; ++i) samples.z[i + j * n] = z[i];
    }
    samples.psi[j] = psi;
    for (std::size_t k = 0; k < n_theta; ++k) {
      samples.theta[k + j * n_theta] = theta[k];
    }

    // Accumulate allocation counts
    for (std::size_t i = 0; i < n; ++i) {
      samples.alloc[i * nc + static_cast<std::size_t>(z[i])] += 1;
    }

    // Sufficient statistic
    samples.sufficient_stat[j] = model.SufficientStatistic(z);

    // Joint MAP tracking
    double score = model.SpatialLogLikelihood(z, psi) + LogPriorPsi(psi);
    if (model.has_emission()) {
      score += EmissionLogLikelihood(y, z, theta, model.emission_type());
    }
    if (score > samples.log_posterior_map) {
      samples.log_posterior_map = score;
      samples.z_map.assign(z.begin(), z.end());
      samples.map_iteration = j;
    }
  }

  return samples;
}

}  // namespace mdgm
