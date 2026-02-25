#pragma once

#include <cstddef>
#include <mdgm/emission.hpp>
#include <mdgm/observations.hpp>
#include <mdgm/rng.hpp>
#include <mdgm/spatial_random_field.hpp>
#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace mdgm {

class Model {
 public:
  // Standalone (direct data model, no emission layer)
  explicit Model(std::unique_ptr<SpatialRandomField> spatial);

  // Hierarchical (hidden field + emission)
  Model(std::unique_ptr<SpatialRandomField> spatial, FamilyType emission);

  bool has_emission() const;
  FamilyType emission_type() const;

  // Full conditional for z_i: spatial * emission (if present), normalized
  std::vector<double> ZFullConditional(
      std::span<const int> z, std::size_t vertex, double psi,
      const Observations& y, std::span<const double> eta) const;

  // Spatial log-likelihood (for beta MH ratio)
  double SpatialLogLikelihood(std::span<const int> z, double psi) const;

  // Update graph (delegates to spatial field)
  void UpdateGraph(std::span<const int> z, double psi, RNG& rng);

  // Update emission params (only valid if has_emission())
  std::vector<double> UpdateEmissionParams(
      const Observations& y, std::span<const int> z,
      std::span<const double> eta, std::span<const double> prior_params,
      RNG& rng) const;

  // Store current graph sample
  void StoreDagSample(std::vector<std::size_t>& dag_data,
                      std::size_t iteration) const;

  std::size_t nvertices() const;
  std::size_t ncolors() const;
  const SpatialRandomField& spatial() const;

 private:
  std::unique_ptr<SpatialRandomField> spatial_;
  std::optional<FamilyType> emission_type_;
};

}  // namespace mdgm
