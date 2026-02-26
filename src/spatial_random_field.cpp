#include <cmath>
#include <cstddef>
#include <mdgm/mcmc.hpp>
#include <mdgm/rng.hpp>
#include <mdgm/spatial_random_field.hpp>
#include <span>

namespace mdgm {

double SpatialRandomField::UpdatePsi(
    std::span<const int> z, double psi,
    double psi_tune, std::size_t& accepted, RNG& rng) {
  double proposal = psi + rng.normal(0.0, psi_tune);
  double log_proposed = LogLikelihood(z, proposal) + LogPriorPsi(proposal);
  double log_current = LogLikelihood(z, psi) + LogPriorPsi(psi);
  if (std::log(rng.uniform()) < log_proposed - log_current) {
    ++accepted;
    return proposal;
  }
  return psi;
}

}  // namespace mdgm
