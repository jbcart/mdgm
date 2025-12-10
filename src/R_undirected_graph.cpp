#include <cpp11.hpp>
#include <cstddef>
#include <mdgm/graph_storage.hpp>
#include <mdgm/rng.hpp>
#include <mdgm/undirected_graph.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

inline mdgm::SpanningTreeMethod ParseMethod(const std::string& method) {
  using namespace mdgm;
  if (method == "wilson") return kWilson;
  if (method == "aldous_broder") return kAldousBroder;
  if (method == "fast_forward") return kFastForward;
  throw std::invalid_argument("Unknown spanning tree method: " + method);
}

inline void CheckLengths(const cpp11::integers& row, const cpp11::integers& col) {
  if (row.size() != col.size()) {
    throw std::invalid_argument("row and col must have the same length");
  }
}

}  // namespace

[[cpp11::register]]
cpp11::external_pointer<mdgm::UndirectedGraph> undirected_graph_create_cpp(
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
  auto graph = std::make_unique<mdgm::UndirectedGraph>(coo);
  return cpp11::external_pointer<mdgm::UndirectedGraph>(graph.release());
}

[[cpp11::register]]
cpp11::writable::list undirected_graph_sample_spanning_tree_cpp(
    cpp11::external_pointer<mdgm::UndirectedGraph> g, cpp11::strings method, int k,
    cpp11::external_pointer<mdgm::RNG> rng) {
  using namespace cpp11::literals;
  auto meth = ParseMethod(static_cast<std::string>(method[0]));

  mdgm::GraphCOO tree = g->SampleSpanningTree(*rng, meth, k);

  const auto& row_ind = tree.row_ind();
  const auto& col_ind = tree.col_ind();

  cpp11::writable::integers r(row_ind.size());
  cpp11::writable::integers c(col_ind.size());
  for (std::size_t i = 0; i < row_ind.size(); ++i) {
    r[i] = static_cast<int>(row_ind[i] + 1);
    c[i] = static_cast<int>(col_ind[i] + 1);
  }
  return cpp11::writable::list({"row"_nm = r, "col"_nm = c});
}

[[cpp11::register]]
cpp11::writable::list undirected_graph_neighbors_cpp(
    cpp11::external_pointer<mdgm::UndirectedGraph> g, int vertex) {
  using namespace cpp11::literals;
  if (vertex < 1 || vertex > static_cast<int>(g->nvertices())) {
    cpp11::stop("vertex must be between 1 and nvertices");
  }
  std::span<const std::size_t> nbrs = g->neighbors(static_cast<std::size_t>(vertex - 1));
  std::span<const double> weights = g->neighbor_weights(static_cast<std::size_t>(vertex - 1));

  cpp11::writable::integers c(nbrs.size());
  cpp11::writable::doubles w(weights.size());
  for (std::size_t i = 0; i < nbrs.size(); ++i) {
    c[i] = static_cast<int>(nbrs[i] + 1);
    w[i] = weights[i];
  }
  return cpp11::writable::list({"neighbors"_nm = c, "weights"_nm = w});
}

[[cpp11::register]]
int undirected_graph_nvertices_cpp(cpp11::external_pointer<mdgm::UndirectedGraph> g) {
  return static_cast<int>(g->nvertices());
}

[[cpp11::register]]
int undirected_graph_nedges_cpp(cpp11::external_pointer<mdgm::UndirectedGraph> g) {
  return static_cast<int>(g->nedges());
}
