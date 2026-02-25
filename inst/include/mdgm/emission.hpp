#pragma once

#include <cstddef>
#include <mdgm/observations.hpp>
#include <mdgm/rng.hpp>
#include <mdgm/spatial_random_field.hpp>
#include <span>
#include <vector>

namespace mdgm {

// Returns likelihood p(y_i | z_i = k, theta) for each color k
std::vector<double> EmissionLikelihood(
    std::span<const double> y_i, std::span<const double> theta,
    std::size_t ncolors, FamilyType type);

// Log-likelihood of all observations: sum_i log p(y_i | z_i, theta)
double EmissionLogLikelihood(
    const Observations& y, std::span<const int> z,
    std::span<const double> theta, FamilyType type);

// Conjugate update for emission parameters
// prior_params interpretation depends on FamilyType:
//   kBernoulli: {a, b} for Beta(a,b) prior on each p_k
std::vector<double> UpdateEmissionParams(
    const Observations& y, std::span<const int> z,
    std::span<const double> theta, std::span<const double> prior_params,
    std::size_t ncolors, FamilyType type, RNG& rng);

}  // namespace mdgm
