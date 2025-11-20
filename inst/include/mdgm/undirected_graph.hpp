#pragma once

#include <cstddef>
#include <mdgm/graph_storage.hpp>
#include <mdgm/rng.hpp>
#include <span>
#include <vector>

namespace mdgm {

enum SpanningTreeMethod {
  kWilson,
  kAldousBroder,
  kHybrid,
  kFastForward,
};

class UndirectedGraph {
 public:
  UndirectedGraph(const GraphCSR& csr);
  UndirectedGraph(const GraphCOO& coo);

  std::span<const std::size_t> neighbors(std::size_t vertex) const;
  std::span<const double> neighbor_weights(std::size_t vertex) const;
  std::size_t nvertices() const noexcept;
  std::size_t nedges() const noexcept;

  GraphCOO SampleSpanningTree(RNG& rng, SpanningTreeMethod method = kWilson, int k = 1000) const;

 private:
  void ValidateUndirected_() const;
  void ValidateConnected_() const;

  GraphCOO SampleSpanningTreeWilson_(RNG& rng) const;
  GraphCOO SampleSpanningTreeAldousBroder_(RNG& rng) const;
  GraphCOO SampleSpanningTreeFastForward_(RNG& rng, int k) const;
  GraphCOO SampleSpanningTreeHybrid_(RNG& rng, int k) const;

  std::size_t NextVertex_(RNG& rng, std::size_t current) const;
  std::size_t SampleRootVertex_(RNG& rng) const;

  GraphCSR csr_;
  [[maybe_unused]] bool is_connected_;
  [[maybe_unused]] bool checked_connected_;
};

}  // namespace mdgm
