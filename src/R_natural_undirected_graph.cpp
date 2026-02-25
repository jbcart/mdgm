#include <cpp11.hpp>
#include <cstddef>
#include <mdgm/directed_acyclic_graph.hpp>
#include <mdgm/graph_storage.hpp>
#include <mdgm/natural_undirected_graph.hpp>
#include <mdgm/rng.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

inline mdgm::SpanningTreeMethod ParseMethod(const std::string& method) {
  using namespace mdgm;
  if (method == "wilson") return SpanningTreeMethod::kWilson;
  if (method == "aldous_broder") return SpanningTreeMethod::kAldousBroder;
  if (method == "fast_forward") return SpanningTreeMethod::kFastForward;
  throw std::invalid_argument("Unknown spanning tree method: " + method);
}

inline void CheckLengths(const cpp11::integers& row, const cpp11::integers& col) {
  if (row.size() != col.size()) {
    throw std::invalid_argument("row and col must have the same length");
  }
}

}  // namespace

[[cpp11::register]]
cpp11::external_pointer<mdgm::NaturalUndirectedGraph> nug_create_cpp(
    int nvertices, const cpp11::integers& row, const cpp11::integers& col,
    const cpp11::doubles& weights) {
  CheckLengths(row, col);
  std::vector<std::size_t> row_ind(row.size());
  std::vector<std::size_t> col_ind(col.size());
  std::vector<double> edge_weights(weights.size());
  for (std::size_t i{0}; i < row.size(); ++i) {
    if (row[i] < 1 || row[i] > nvertices || col[i] < 1 || col[i] > nvertices) {
      cpp11::stop("row and col values must be between 1 and nvertices");
    }
    row_ind[i] = static_cast<std::size_t>(row[i] - 1);
    col_ind[i] = static_cast<std::size_t>(col[i] - 1);
    edge_weights[i] = weights[i];
  }
  mdgm::GraphCOO coo(nvertices, row_ind, col_ind, edge_weights);
  auto graph = std::make_unique<mdgm::NaturalUndirectedGraph>(coo);
  return cpp11::external_pointer<mdgm::NaturalUndirectedGraph>(graph.release());
}

[[cpp11::register]]
cpp11::writable::list nug_sample_spanning_tree_cpp(
    cpp11::external_pointer<mdgm::NaturalUndirectedGraph> g, cpp11::strings method, int k,
    cpp11::external_pointer<mdgm::RNG> rng) {
  using namespace cpp11::literals;
  auto meth = ParseMethod(static_cast<std::string>(method[0]));

  mdgm::DirectedAcyclicGraph dag = g->SampleSpanningTree(*rng, meth, k);

  // Build edge list from DAG parent structure
  std::vector<int> parents_r;
  std::vector<int> children_r;
  for (std::size_t i = 0; i < dag.nvertices(); ++i) {
    for (std::size_t p : dag.parents(i)) {
      parents_r.push_back(static_cast<int>(p + 1));
      children_r.push_back(static_cast<int>(i + 1));
    }
  }

  cpp11::writable::integers r(parents_r.size());
  cpp11::writable::integers c(children_r.size());
  for (std::size_t i = 0; i < parents_r.size(); ++i) {
    r[i] = parents_r[i];
    c[i] = children_r[i];
  }
  return cpp11::writable::list({"parent"_nm = r, "child"_nm = c});
}

[[cpp11::register]]
cpp11::writable::integers nug_neighbors_cpp(
    cpp11::external_pointer<mdgm::NaturalUndirectedGraph> g, int vertex) {
  if (vertex < 1 || vertex > static_cast<int>(g->nvertices())) {
    cpp11::stop("vertex must be between 1 and nvertices");
  }
  std::span<const std::size_t> nbrs = g->neighbors(static_cast<std::size_t>(vertex - 1));

  cpp11::writable::integers c(nbrs.size());
  for (std::size_t i = 0; i < nbrs.size(); ++i) {
    c[i] = static_cast<int>(nbrs[i] + 1);
  }
  return c;
}

[[cpp11::register]]
cpp11::writable::doubles nug_neighbor_weights_cpp(
    cpp11::external_pointer<mdgm::NaturalUndirectedGraph> g, int vertex) {
  if (vertex < 1 || vertex > static_cast<int>(g->nvertices())) {
    cpp11::stop("vertex must be between 1 and nvertices");
  }
  std::span<const double> weights = g->neighbor_weights(static_cast<std::size_t>(vertex - 1));

  cpp11::writable::doubles w(weights.size());
  for (std::size_t i = 0; i < weights.size(); ++i) {
    w[i] = weights[i];
  }
  return w;
}

[[cpp11::register]]
cpp11::external_pointer<mdgm::NaturalUndirectedGraph> nug_generate_regular_cpp(
    int nrows, int ncols, int order) {
  std::vector<std::size_t> dims = {static_cast<std::size_t>(nrows),
                                    static_cast<std::size_t>(ncols)};
  auto graph = std::make_unique<mdgm::NaturalUndirectedGraph>(
      mdgm::GenerateRegularGraph(dims, order));
  return cpp11::external_pointer<mdgm::NaturalUndirectedGraph>(graph.release());
}

[[cpp11::register]]
int nug_nvertices_cpp(cpp11::external_pointer<mdgm::NaturalUndirectedGraph> g) {
  return static_cast<int>(g->nvertices());
}

[[cpp11::register]]
int nug_nedges_cpp(cpp11::external_pointer<mdgm::NaturalUndirectedGraph> g) {
  return static_cast<int>(g->nedges());
}
