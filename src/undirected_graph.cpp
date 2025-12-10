#include <algorithm>
#include <mdgm/graph_storage.hpp>
#include <mdgm/rng.hpp>
#include <mdgm/undirected_graph.hpp>
#include <stdexcept>
#include <vector>

namespace mdgm {

UndirectedGraph GenerateRegularGraph(std::vector<std::size_t> dims, int order) {
  if (dims.size() != 2) {
    throw std::invalid_argument("Only 2D regular graphs are supported");
  }
  if (order < 1 || order > 2) {
    throw std::invalid_argument("Order must be between 1 and 2");
  }
  std::size_t nrows{dims[0]};
  std::size_t ncols{dims[1]};
  std::size_t nvertices{nrows * ncols};
  std::size_t nedges{};
  if (order == 1) {
    nedges = nvertices * 2 - nrows - ncols;
  } else if (order == 2) {
    nedges = nvertices * 3 - nrows * 2 - ncols * 2 + 1;
  }
  std::vector<std::size_t> row_ptr(nvertices + 1, 0);
  std::vector<std::size_t> col_ind;
  col_ind.reserve(nedges);
  std::vector<double> weights(nedges, 1.0);
  std::size_t v{};
  for (std::size_t r{0}; r < nrows; ++r) {
    for (std::size_t c{0}; c < ncols; ++c) {
      v = r * ncols + c;
      row_ptr[v + 1] = row_ptr[v];
      // add SW neighbor
      if (order == 2 && r > 0 && c > 0) {
        col_ind.push_back(v - ncols - 1);
        ++row_ptr[v + 1];
      }
      // add S neighbor
      if (r > 0) {
        col_ind.push_back(v - ncols);
        ++row_ptr[v + 1];
      }
      // add SE neighbor
      if (order == 2 && r > 0 && c + 1 < ncols) {
        col_ind.push_back(v - ncols + 1);
        ++row_ptr[v + 1];
      }
      // add W neighbor
      if (c > 0) {
        col_ind.push_back(v - 1);
        ++row_ptr[v + 1];
      }
      // add E neighbor
      if (c + 1 < ncols) {
        col_ind.push_back(v + 1);
        ++row_ptr[v + 1];
      }
      // add NW neighbor
      if (order == 2 && r + 1 < nrows && c > 0) {
        col_ind.push_back(v + ncols - 1);
        ++row_ptr[v + 1];
      }
      // add N neighbor
      if (r + 1 < nrows) {
        col_ind.push_back(v + ncols);
        ++row_ptr[v + 1];
      }
      // add NE neighbor
      if (order == 2 && r + 1 < nrows && c + 1 < ncols) {
        col_ind.push_back(v + ncols + 1);
        ++row_ptr[v + 1];
      }
    }
  }
  return UndirectedGraph(GraphCSR(nvertices, row_ptr, col_ind, weights), false);
}

UndirectedGraph::UndirectedGraph(const GraphCSR& csr, bool validate) : csr_(csr) {
  if (validate) {
    ValidateUndirected_();
    ValidateConnected_();
  }
}

UndirectedGraph::UndirectedGraph(const GraphCOO& coo) : csr_(coo) {
  ValidateUndirected_();
  ValidateConnected_();
}

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
  if (!is_connected_) {
    throw std::runtime_error("Wilson's algorithm requires connected graph");
  }
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
  if (!is_connected_) {
    throw std::runtime_error("Aldous-Broder algorithm requires connected graph");
  }
  std::vector<std::size_t> row_ind;
  std::vector<std::size_t> col_ind;
  const std::size_t n = nvertices();
  row_ind.reserve(n - 1);
  col_ind.reserve(n - 1);
  std::vector<uint8_t> in_tree(n, 0);
  std::size_t in_tree_count{0};

  // initialize by sampling root vertex from stationary distribution
  std::size_t current = SampleRootVertex_(rng);
  in_tree[current] = 1;
  ++in_tree_count;

  // perform random walk until all vertices are in the tree
  while (in_tree_count < n) {
    std::size_t next = NextVertex_(rng, current);
    if (!in_tree[next]) {
      // add edge to tree
      row_ind.push_back(current);
      col_ind.push_back(next);
      in_tree[next] = 1;
      ++in_tree_count;
    }
    current = next;
  }
  if (row_ind.size() != n - 1) {
    throw std::logic_error("Internal error in Aldous-Broder algorithm: incorrect number of edges");
  }
  return GraphCOO(nvertices(), row_ind, col_ind);
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
          throw std::invalid_argument("Graph is not undirected: missing reverse edge (row_ind = " +
                                      std::to_string(v) + ", col_ind = " + std::to_string(u) + ")");
        }
      }
    }
  }
}

void UndirectedGraph::ValidateConnected_() {
  const std::size_t n = csr_.nvertices();
  std::vector<bool> visited(n, false);
  std::vector<std::size_t> stack;
  stack.push_back(0);
  visited[0] = true;
  std::size_t visit_count{1};

  while (!stack.empty()) {
    std::size_t current = stack.back();
    stack.pop_back();
    std::span<const std::size_t> nbrs = csr_.adjacent(current);
    for (std::size_t nbr : nbrs) {
      if (!visited[nbr]) {
        visited[nbr] = true;
        ++visit_count;
        stack.push_back(nbr);
      }
    }
  }

  if (visit_count < n) {
    is_connected_ = false;
    std::fprintf(stderr, "Warning: Graph is not connected; spanning tree algorithms may fail.\n");
  } else {
    is_connected_ = true;
  }
}

}  // namespace mdgm
