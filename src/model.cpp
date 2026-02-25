#include <cstddef>
#include <mdgm/emission.hpp>
#include <mdgm/model.hpp>
#include <mdgm/observations.hpp>
#include <mdgm/rng.hpp>
#include <mdgm/spatial_random_field.hpp>
#include <span>
#include <stdexcept>
#include <vector>

namespace mdgm {

Model::Model(std::unique_ptr<SpatialRandomField> spatial)
    : spatial_(std::move(spatial)), emission_type_(std::nullopt) {}

Model::Model(std::unique_ptr<SpatialRandomField> spatial, FamilyType emission)
    : spatial_(std::move(spatial)), emission_type_(emission) {}

bool Model::has_emission() const { return emission_type_.has_value(); }

FamilyType Model::emission_type() const {
  if (!emission_type_.has_value()) {
    throw std::logic_error("No emission type on standalone model");
  }
  return *emission_type_;
}

std::vector<double> Model::ZFullConditional(
    std::span<const int> z, std::size_t vertex, double psi,
    const Observations& y, std::span<const double> theta) const {
  auto probs = spatial_->ZFullConditional(z, vertex, psi);

  if (has_emission() && !y.empty(vertex)) {
    auto lik = EmissionLikelihood(y[vertex], theta, spatial_->ncolors(),
                                  *emission_type_);
    for (std::size_t k = 0; k < probs.size(); ++k) {
      probs[k] *= lik[k];
    }
    // Renormalize
    double sum = 0.0;
    for (double p : probs) sum += p;
    if (sum > 0.0) {
      for (double& p : probs) p /= sum;
    }
  }

  return probs;
}

double Model::SpatialLogLikelihood(std::span<const int> z, double psi) const {
  return spatial_->LogLikelihood(z, psi);
}

void Model::UpdateGraph(std::span<const int> z, double psi, RNG& rng) {
  spatial_->UpdateGraph(z, psi, rng);
}

std::vector<double> Model::UpdateEmissionParams(
    const Observations& y, std::span<const int> z,
    std::span<const double> theta, std::span<const double> prior_params,
    RNG& rng) const {
  if (!has_emission()) {
    throw std::logic_error("Cannot update emission params on standalone model");
  }
  return mdgm::UpdateEmissionParams(y, z, theta, prior_params,
                                     spatial_->ncolors(), *emission_type_, rng);
}

void Model::StoreDagSample(std::vector<std::size_t>& dag_data,
                           std::size_t iteration) const {
  spatial_->StoreSample(dag_data, iteration, spatial_->nvertices());
}

std::size_t Model::nvertices() const { return spatial_->nvertices(); }
std::size_t Model::ncolors() const { return spatial_->ncolors(); }
const SpatialRandomField& Model::spatial() const { return *spatial_; }

}  // namespace mdgm
