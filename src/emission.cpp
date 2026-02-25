#include <cmath>
#include <cstddef>
#include <mdgm/emission.hpp>
#include <mdgm/observations.hpp>
#include <mdgm/rng.hpp>
#include <mdgm/spatial_random_field.hpp>
#include <span>
#include <stdexcept>
#include <vector>

namespace mdgm {

std::vector<double> EmissionLikelihood(
    std::span<const int> y_i, std::span<const double> eta,
    std::size_t ncolors, FamilyType type) {
  std::vector<double> lik(ncolors, 1.0);

  switch (type) {
    case FamilyType::kBernoulli:
      for (std::size_t k = 0; k < ncolors; ++k) {
        double p = eta[k];
        for (int y : y_i) {
          lik[k] *= (y == 1) ? p : (1.0 - p);
        }
      }
      break;
    default:
      throw std::invalid_argument("Unsupported emission family type");
  }

  return lik;
}

double EmissionLogLikelihood(
    const Observations& y, std::span<const int> z,
    std::span<const double> eta, FamilyType type) {
  double ll = 0.0;

  switch (type) {
    case FamilyType::kBernoulli:
      for (std::size_t i = 0; i < y.nvertices(); ++i) {
        if (y.empty(i)) continue;
        double p = eta[static_cast<std::size_t>(z[i])];
        for (int yij : y[i]) {
          ll += (yij == 1) ? std::log(p) : std::log(1.0 - p);
        }
      }
      break;
    default:
      throw std::invalid_argument("Unsupported emission family type");
  }

  return ll;
}

std::vector<double> UpdateEmissionParams(
    const Observations& y, std::span<const int> z,
    std::span<const double> eta, std::span<const double> prior_params,
    std::size_t ncolors, FamilyType type, RNG& rng) {
  std::vector<double> new_eta(eta.begin(), eta.end());

  switch (type) {
    case FamilyType::kBernoulli: {
      double a = prior_params[0];
      double b = prior_params[1];

      std::vector<double> successes(ncolors, 0.0);
      std::vector<double> trials(ncolors, 0.0);
      for (std::size_t i = 0; i < y.nvertices(); ++i) {
        std::size_t k = static_cast<std::size_t>(z[i]);
        for (int yij : y[i]) {
          trials[k] += 1.0;
          successes[k] += static_cast<double>(yij);
        }
      }

      // Truncated Beta posteriors with identifiability: eta[0] < eta[1]
      new_eta[0] = rng.truncated_beta(
          a + successes[0], b + trials[0] - successes[0],
          0.0, new_eta[1]);
      new_eta[1] = rng.truncated_beta(
          a + successes[1], b + trials[1] - successes[1],
          new_eta[0], 1.0);
      break;
    }
    default:
      throw std::invalid_argument("Unsupported emission family type");
  }

  return new_eta;
}

}  // namespace mdgm
