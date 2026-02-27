#pragma once

#include <cstddef>
#include <mdgm/rng.hpp>
#include <span>
#include <vector>

namespace mdgm {

enum class FamilyType {
  kCategorical,
  kGaussian,
  kPoisson,
  kBernoulli,
};

class SpatialRandomField {
 public:
  virtual ~SpatialRandomField() = default;

  // Spatial full conditional: p(z_i = k | z_{-i}, psi) for each color k
  // Pure spatial component — does NOT include emission
  virtual std::vector<double> ZFullConditional(
      std::span<const int> z, std::size_t vertex, double psi) const = 0;

  // Log-likelihood of z under this spatial model (pseudo-log-likelihood)
  virtual double LogLikelihood(std::span<const int> z, double psi) const = 0;

  // Update graph structure given current z and psi (e.g., sample new DAG)
  virtual void UpdateGraph(std::span<const int> z, double psi, RNG& rng) = 0;

  // Store compact graph representation for posterior storage
  virtual void StoreSample(std::vector<std::size_t>& dag_data,
                           std::size_t iteration, std::size_t n) const = 0;

  // Update psi via MH step. Default: random walk with LogLikelihood + Half-Cauchy prior.
  // Subclasses (e.g., MRF exchange algorithm) may override.
  virtual double UpdatePsi(std::span<const int> z, double psi,
                           double psi_tune, std::size_t& accepted, RNG& rng);

  // Sufficient statistic: count of same-color neighbor pairs (each edge once)
  virtual double SufficientStatistic(std::span<const int> z) const = 0;

  virtual std::size_t nvertices() const = 0;
  virtual std::size_t ncolors() const = 0;
};

}  // namespace mdgm
