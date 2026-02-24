#pragma once

#include <cstddef>
#include <mdgm/graph_storage.hpp>
#include <mdgm/directed_acyclic_graph.hpp>
#include <mdgm/rng.hpp>
#include <span>
#include <vector>

namespace mdgm {

enum SpanningTreeMethod {
  kWilson,
  kAldousBroder,
  kFastForward,
};

class NaturalUndirectedGraph {
 public:
  NaturalUndirectedGraph(const GraphCOO& coo);
  NaturalUndirectedGraph(const GraphCSR& csr, bool validate = true);

  ~NaturalUndirectedGraph() = default;
  NaturalUndirectedGraph(const NaturalUndirectedGraph&) = default;
  NaturalUndirectedGraph& operator=(const NaturalUndirectedGraph&) = default;

  std::span<const std::size_t> neighbors(std::size_t vertex) const;
  std::span<const double> neighbor_weights(std::size_t vertex) const;
  std::size_t nvertices() const noexcept;
  std::size_t nedges() const noexcept;

  DirectedAcyclicGraph SampleSpanningTree(RNG& rng, SpanningTreeMethod method = kWilson, int k = 1000) const;
  DirectedAcyclicGraph SampleAcyclicOrientation(RNG& rng) const;

 private:
  void ValidateUndirected_() const;
  void ValidateConnected_();

  DirectedAcyclicGraph SampleSpanningTreeWilson_(RNG& rng) const;
  DirectedAcyclicGraph SampleSpanningTreeAldousBroder_(RNG& rng) const;
  DirectedAcyclicGraph SampleSpanningTreeFastForward_(RNG& rng, int k) const;

  std::size_t NextVertex_(RNG& rng, std::size_t current) const;
  std::size_t SampleRootVertex_(RNG& rng) const;
  void UpdateRemaining_(std::vector<std::size_t>& remaining,
                        std::vector<std::size_t>& pos, std::size_t vertex) const;

  GraphCSR csr_;
  bool is_connected_;
};

NaturalUndirectedGraph GenerateRegularGraph(std::vector<std::size_t> dims, int order);

}  // namespace mdgm
