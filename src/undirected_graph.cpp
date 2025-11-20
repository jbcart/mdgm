#include <algorithm>
#include <mdgm/graph_storage.hpp>
#include <mdgm/rng.hpp>
#include <mdgm/undirected_graph.hpp>
#include <stdexcept>
#include <vector>

namespace mdgm {

UndirectedGraph::UndirectedGraph(const GraphCSR& csr) : csr_(csr) { ValidateUndirected_(); }

UndirectedGraph::UndirectedGraph(const GraphCOO& coo) : csr_(coo) { ValidateUndirected_(); }

std::span<const std::size_t> UndirectedGraph::neighbors(std::size_t vertex) const {
  return csr_.adjacent(vertex);
}

std::span<const double> UndirectedGraph::neighbor_weights(std::size_t vertex) const {
  return csr_.adjacent_weights(vertex);
}

std::size_t UndirectedGraph::nvertices() const noexcept { return csr_.nvertices(); }

std::size_t UndirectedGraph::nedges() const noexcept { return csr_.weights().size() / 2; }

GraphCOO UndirectedGraph::SampleSpanningTree(RNG& rng, SpanningTreeMethod method, int k) const {
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
  // NOTE: This samples a root vertex from the stationary distribution of the random walk along the
  // edges of the graph.
  return csr_.col_ind()[rng.discrete(csr_.weights())];
}

void UndirectedGraph::UpdateRemaining_(std::vector<std::size_t>& remaining,
                                       std::vector<std::size_t>& pos, std::size_t vertex) const {
  // swap with last and pop, update pos
  if (pos[vertex] == SIZE_MAX) return;  // already removed
  std::size_t idx = pos[vertex];
  std::size_t last = remaining.back();
  remaining[idx] = last;
  pos[last] = idx;
  remaining.pop_back();
  pos[vertex] = SIZE_MAX;
}

GraphCOO UndirectedGraph::SampleSpanningTreeWilson_(RNG& rng) const {
  // store as COO representation
  std::vector<std::size_t> row_ind;
  std::vector<std::size_t> col_ind;
  const std::size_t n = nvertices();
  row_ind.reserve(n - 1);
  col_ind.reserve(n - 1);

  std::vector<int> in_tree(n, 0);
  std::size_t in_tree_count{0};

  // current random walk with augmented data for loop erasure
  std::vector<std::size_t> walk;
  walk.reserve(n);
  std::vector<std::size_t> visit_time(n, SIZE_MAX);

  // keeping track of vertices not yet in the tree for efficient sampling
  std::vector<std::size_t> remaining(n);
  std::vector<std::size_t> pos(n, SIZE_MAX);
  for (std::size_t i = 0; i < n; ++i) {
    remaining[i] = i;
    pos[i] = i;
  }

  // initialize by sampling root vertex from stationary distribution
  std::size_t root = SampleRootVertex_(rng);
  UpdateRemaining_(remaining, pos, root);

  in_tree[root] = 1;
  ++in_tree_count;
  
  // iterate while there exists at least one vertex not yet in the tree
  while (in_tree_count < n) {
    walk.clear();
    // start a new random walk from a random vertex not in the tree
    std::size_t start = remaining[rng.uniform<std::size_t>(0, remaining.size() - 1)];
    walk.push_back(start);
    visit_time[start] = 0;

    std::size_t next = NextVertex_(rng, start);
    while (in_tree[next] == 0) {
      if (visit_time[next] != SIZE_MAX) {
        // loop detected, erase loop
        std::size_t loop_start = visit_time[next];
        while (walk.size() > loop_start + 1) {
          visit_time[walk.back()] = SIZE_MAX;
          walk.pop_back();
        }
      } else {
        visit_time[next] = walk.size();
        walk.push_back(next);
      }
      next = NextVertex_(rng, next);
    }

    // add final vertex that connects to tree
    walk.push_back(next);

    // add the walk to the tree (edges directed away from root), mark vertices in tree
    for (std::size_t i = 0; i < walk.size() - 1; ++i) {
      std::size_t child = walk[i];
      std::size_t parent = walk[i + 1];
      if (!in_tree[child]) {
        UpdateRemaining_(remaining, pos, child);
        in_tree[child] = 1;
        ++in_tree_count;
      }
      visit_time[child] = SIZE_MAX;  // reset visit time
      row_ind.push_back(parent);
      col_ind.push_back(child);
    }
  }
  if (row_ind.size() != n - 1) {
    throw std::logic_error("Internal error in Wilson's algorithm: incorrect number of edges");
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

GraphCOO UndirectedGraph::SampleSpanningTreeFastForward_(RNG& rng, int k) const {
  // unimplemented
  return GraphCOO(0, {}, {}, {});
}

void UndirectedGraph::ValidateUndirected_() const {
  const std::size_t n = csr_.nvertices();
  for (std::size_t u = 0; u < n; ++u) {
    std::span<const std::size_t> nbrs_u = csr_.adjacent(u);
    for (std::size_t v : nbrs_u) {
      if (v == u) {
        throw std::invalid_argument("Self-loops are not allowed in undirected graph");
      }
      if (u < v) {
        std::span<const std::size_t> nbrs_v = csr_.adjacent(v);
        auto it = std::lower_bound(nbrs_v.begin(), nbrs_v.end(), u);
        if (it == nbrs_v.end() || *it != u) {
          throw std::invalid_argument(
              "Graph is not undirected: missing reverse edge (row_ind = " + std::to_string(v) +
              ", col_ind = " + std::to_string(u) + ")");
        }
      }
    }
  }
}

void UndirectedGraph::ValidateConnected_() const {
  // unimplemented
}

}  // namespace mdgm
