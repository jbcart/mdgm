#include "undirected_graph.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "graph_storage.h"
#include "random.h"

namespace mdgm {

UndirectedGraph::UndirectedGraph(const GraphCSR& csr) : csr_(csr) {
  ValidateUndirected_();
}

UndirectedGraph::UndirectedGraph(const GraphCOO& coo) : csr_(coo) {
  ValidateUndirected_();
}

std::span<const std::size_t> UndirectedGraph::neighbors(
    std::size_t vertex) const {
  return csr_.adjacent(vertex);
}

std::span<const double> UndirectedGraph::neighbor_weights(
    std::size_t vertex) const {
  return csr_.adjacent_weights(vertex);
}

std::size_t UndirectedGraph::nvertices() const noexcept {
  return csr_.nvertices();
}

std::size_t UndirectedGraph::nedges() const noexcept {
  return csr_.weights().size() / 2;
}

GraphCOO UndirectedGraph::SampleSpanningTree(
    RNG& rng, SpanningTreeMethod method = kWilson, int k = 1000) const {
  switch (method) {
    case kWilson:
      return SampleSpanningTreeWilson_(rng);
    case kAldousBroder:
      return SampleSpanningTreeAldousBroder_(rng);
    case kHybrid:
      return SampleSpanningTreeHybrid_(rng, k);
    case kFastForward:
      return SampleSpanningTreeFastForward_(rng, k);
    default:
      throw std::invalid_argument("Unknown spanning tree method");
  }
}

std::size_t UndirectedGraph::NextVertex_(RNG& rng, std::size_t current) const {
  std::span<const double> nbr_weights = neighbor_weights(current);
  std::span<const std::size_t> nbrs = neighbors(current);
  return nbrs[rng.discrete(nbr_weights)];
}

std::size_t UndirectedGraph::SampleRootVertex_(RNG& rng) const {
  return csr_.col_ind()[rng.discrete(csr_.weights())];
}

GraphCOO UndirectedGraph::SampleSpanningTreeWilson_(RNG& rng) const {
  // store as COO representation
  std::vector<std::size_t> row_ind;
  std::vector<std::size_t> col_ind;
  row_ind.reserve(nvertices() - 1);
  col_ind.reserve(nvertices() - 1);

  std::vector<int> in_tree(nvertices(), 0);
  std::vector<std::size_t> random_walk;
  std::size_t root = SampleRootVertex_(rng);

  in_tree[root] = 1;
  random_walk.reserve(nvertices());

  // iterate while there exists at least one vertex not yet in the tree
  while (std::count(in_tree.begin(), in_tree.end(), 1) < nvertices()) {
    // start a new random walk from a random vertex not in the tree
    std::size_t next = rng.uniform<std::size_t>(
        0, std::count(in_tree.begin(), in_tree.end(), 0) - 1);
    random_walk.push_back(next +
                          std::count(in_tree.data(), in_tree.data() + next, 1));

    next = NextVertex_(rng, next);
    while (in_tree[next] == 0) {
      auto loop_start = std::find(random_walk.begin(), random_walk.end(), next);
      if (loop_start != random_walk.end()) {
        random_walk.erase(loop_start, random_walk.end());
      }
      random_walk.push_back(next);
      next = NextVertex_(rng, next);
    }
    random_walk.push_back(next);
    // add the walk to the tree (edges directed away from root)
    for (std::size_t i = 0; i + 1 < random_walk.size(); ++i) {
      std::size_t to = random_walk[i];
      std::size_t from = random_walk[i + 1];
      row_ind.push_back(to);
      col_ind.push_back(from);
      in_tree[to] = 1;
    }
    random_walk.clear();
  }
  return GraphCOO(nvertices(), row_ind, col_ind);
}

GraphCOO UndirectedGraph::SampleSpanningTreeAldousBroder_(RNG& rng) const {
  // unimplemented
  return GraphCOO(0, {}, {}, {});
}

GraphCOO UndirectedGraph::SampleSpanningTreeHybrid_(RNG& rng, int k) const {
  // unimplemented
  return GraphCOO(0, {}, {}, {});
}

GraphCOO UndirectedGraph::SampleSpanningTreeFastForward_(RNG& rng,
                                                         int k) const {
  // unimplemented
  return GraphCOO(0, {}, {}, {});
}

void UndirectedGraph::ValidateUndirected_() const {
  const std::size_t n = csr_.nvertices();
  for (std::size_t u = 0; u < n; ++u) {
    std::span<const std::size_t> nbrs_u = csr_.adjacent(u);
    for (std::size_t v : nbrs_u) {
      if (v == u) {
        throw std::invalid_argument(
            "Self-loops are not allowed in undirected graph");
      }
      std::span<const std::size_t> nbrs_v = csr_.adjacent(v);
      auto it = std::lower_bound(nbrs_v.begin(), nbrs_v.end(), u);
      if (it == nbrs_v.end()) {
        throw std::invalid_argument(
            "Graph is not undirected: missing reverse edge (" +
            std::to_string(v) + " -> " + std::to_string(u) + ")");
      }
    }
  }
}

void UndirectedGraph::ValidateConnected_() const {
  // unimplemented
}

}  // namespace mdgm
