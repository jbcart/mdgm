#pragma once

#include <cstddef>
#include <mdgm/directed_acyclic_graph.hpp>
#include <mdgm/natural_undirected_graph.hpp>
#include <mdgm/spatial_random_field.hpp>
#include <span>
#include <vector>

namespace mdgm {

class MixtureDirectedGraphicalModels : public SpatialRandomField {
 public:
  MixtureDirectedGraphicalModels(NaturalUndirectedGraph nug, DagType dag_type,
                                  std::size_t n_colors = 2);

  ~MixtureDirectedGraphicalModels() = default;

  // Dispatches on dag_type_:
  //   kSpanningTree: simplified (no v-structures -> neighbors only on tree)
  //   kAcyclicOrientation/kRooted: full conditional with children's normalizers
  std::vector<double> ZFullConditional(
      std::span<const int> z, std::size_t vertex, double psi) const override;

  double LogLikelihood(std::span<const int> z, double psi) const override;

  // kSpanningTree: update weights then sample via Wilson's
  // kAcyclicOrientation: MH with random permutation proposal
  // kRooted: MH with root random walk proposal
  void UpdateGraph(std::span<const int> z, double psi, RNG& rng) override;

  void StoreSample(std::vector<std::size_t>& dag_data,
                   std::size_t iteration, std::size_t n) const override;

  double SufficientStatistic(std::span<const int> z) const override;

  std::size_t nvertices() const override;
  std::size_t ncolors() const override;

  const DirectedAcyclicGraph& current_dag() const;

 private:
  NaturalUndirectedGraph nug_;
  DirectedAcyclicGraph current_dag_;
  DagType dag_type_;
  std::size_t ncolors_;
  std::vector<double> alpha_;  // marginal probs (fixed at 0 for now, future: X*beta)
};

}  // namespace mdgm
