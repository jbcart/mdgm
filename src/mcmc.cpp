#include <cmath>
#include <cstddef>
#include <mdgm/mcmc.hpp>
#include <mdgm/model.hpp>
#include <mdgm/observations.hpp>
#include <mdgm/rng.hpp>
#include <span>
#include <vector>

namespace mdgm {

namespace {

// Half-Cauchy log-prior: log(2 / (pi * (1 + psi^2))) for psi > 0
double LogPriorPsi(double psi) {
  if (psi <= 0.0) return -1e300;
  return std::log(2.0 / M_PI) - std::log(1.0 + psi * psi);
}

}  // namespace

McmcSamples RunMcmc(
    Model& model,
    const Observations& y,
    const std::vector<int>& z_init,
    double psi_init,
    const std::vector<double>& eta_init,
    const McmcConfig& config,
    RNG& rng) {
  const std::size_t n = model.nvertices();
  const std::size_t nc = model.ncolors();
  const std::size_t J = config.n_iterations;

  // Initialize eta
  std::vector<double> eta(eta_init);
  if (eta.empty() && model.has_emission()) {
    eta.resize(nc, 0.5);  // default if not provided
  }

  // Number of emission parameters per iteration
  const std::size_t n_eta = eta.size();

  // Allocate output storage
  McmcSamples samples;
  samples.n_iterations = J;
  samples.n_vertices = n;
  samples.n_colors = nc;
  samples.n_eta = n_eta;
  samples.z.resize(n * J);
  samples.psi.resize(J);
  samples.eta.resize(model.has_emission() ? n_eta * J : 0);
  samples.dag_data.resize(n * J);
  samples.psi_accepted = 0;
  samples.graph_accepted = 0;

  // Initialize iteration 0
  std::vector<int> z(z_init);
  double psi = psi_init;

  for (std::size_t i = 0; i < n; ++i) samples.z[i] = z[i];
  samples.psi[0] = psi;
  for (std::size_t k = 0; k < n_eta; ++k) samples.eta[k] = eta[k];
  model.StoreDagSample(samples.dag_data, 0);

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
        auto probs = model.ZFullConditional(z, i, psi, y, eta);
        z[i] = static_cast<int>(rng.discrete(probs));
      }
    }

    // Step 3: Update psi (MH random walk)
    {
      double proposal = psi + rng.normal(0.0, config.psi_tune);
      double log_proposed = model.SpatialLogLikelihood(z, proposal) +
                            LogPriorPsi(proposal);
      double log_current = model.SpatialLogLikelihood(z, psi) +
                           LogPriorPsi(psi);
      if (std::log(rng.uniform()) < log_proposed - log_current) {
        psi = proposal;
        ++samples.psi_accepted;
      }
    }

    // Step 4: Update emission params (only if hierarchical)
    if (model.has_emission()) {
      eta = model.UpdateEmissionParams(y, z, eta, config.emission_prior_params,
                                        rng);
    }

    // Store samples for this iteration
    for (std::size_t i = 0; i < n; ++i) samples.z[i + j * n] = z[i];
    samples.psi[j] = psi;
    for (std::size_t k = 0; k < n_eta; ++k) {
      samples.eta[k + j * n_eta] = eta[k];
    }
  }

  return samples;
}

}  // namespace mdgm
