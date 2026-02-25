#pragma once

#include <cstddef>
#include <mdgm/graph_storage.hpp>
#include <span>
#include <vector>

namespace mdgm {

enum class DagType {
  kSpanningTree,
  kAcyclicOrientation,
  kRooted,
};

class DirectedAcyclicGraph {
 public:
  DirectedAcyclicGraph(const GraphCOO& coo, DagType type);

  ~DirectedAcyclicGraph() = default;
  DirectedAcyclicGraph(const DirectedAcyclicGraph&) = default;
  DirectedAcyclicGraph& operator=(const DirectedAcyclicGraph&) = default;

  std::span<const std::size_t> children(std::size_t vertex) const;
  std::span<const std::size_t> parents(std::size_t vertex) const;

  std::size_t nvertices() const noexcept;
  std::size_t nedges() const noexcept;

 private:
  GraphCSR csr_;
  GraphCSR transpose_;
  DagType type_;
};

}  // namespace mdgm
