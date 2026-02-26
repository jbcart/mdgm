#pragma once

#include <cstddef>
#include <mdgm/natural_undirected_graph.hpp>
#include <mdgm/rng.hpp>
#include <mdgm/spatial_random_field.hpp>
#include <span>
#include <vector>

namespace mdgm {

enum class PsiMethod {
  kExchange,
  kPseudoLikelihood,
};

class MarkovRandomField : public SpatialRandomField {
 public:
  MarkovRandomField(NaturalUndirectedGraph nug, PsiMethod method,
                    std::size_t n_colors = 2,
                    std::size_t n_aux_sweeps = 200);

  ~MarkovRandomField() = default;

  // Potts full conditional: p(z_i = k | z_{-i}, psi) proportional to
  // exp(psi * count_same_color_neighbors(i, k, z))
  std::vector<double> ZFullConditional(
      std::span<const int> z, std::size_t vertex, double psi) const override;

  // Pseudo-log-likelihood: sum_i log p(z_i | z_{-i}, psi)
  double LogLikelihood(std::span<const int> z, double psi) const override;

  // No-op: MRF has a fixed graph
  void UpdateGraph(std::span<const int> z, double psi, RNG& rng) override;

  // Fills with SIZE_MAX (becomes NA in R)
  void StoreSample(std::vector<std::size_t>& dag_data,
                   std::size_t iteration, std::size_t n) const override;

  // Override for exchange algorithm; pseudo-likelihood uses default
  double UpdatePsi(std::span<const int> z, double psi,
                   double psi_tune, std::size_t& accepted,
                   RNG& rng) override;

  std::size_t nvertices() const override;
  std::size_t ncolors() const override;

  // Draw an exact (binary, CFTP) or approximate (Potts, Gibbs) sample.
  std::vector<int> Sample(double psi, RNG& rng) const;

  // Single CFTP sweep: deterministic update of config using shared uniforms.
  // uniforms is a flat vector of size n * total_sweeps; col_offset selects
  // which sweep's column of uniforms to use.
  void CftpSweep(std::vector<int>& config, double psi,
                 const std::vector<double>& uniforms,
                 std::size_t col_offset) const;

 private:
  // Sufficient statistic: count of same-color neighbor pairs (edge-based)
  double SufficientStatistic(std::span<const int> z) const;

  // Gibbs sample from MRF at given psi (for exchange algorithm auxiliary field)
  std::vector<int> GibbsSample(double psi, RNG& rng) const;

  // CFTP (Propp-Wilson) exact sample from binary MRF via monotone coupling.
  // Only valid for ncolors_ == 2.
  std::vector<int> CftpSample(double psi, RNG& rng) const;

  NaturalUndirectedGraph nug_;
  PsiMethod method_;
  std::size_t ncolors_;
  std::size_t n_aux_sweeps_;
};

}  // namespace mdgm
