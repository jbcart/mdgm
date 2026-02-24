#include <mdgm/directed_acyclic_graph.hpp>
#include <mdgm/graph_storage.hpp>
#include <span>
#include <vector>

namespace mdgm {

DirectedAcyclicGraph::DirectedAcyclicGraph(const GraphCOO& coo, DAGType type)
    : csr_(coo), transpose_(coo.Transpose()), type_(type) {}

std::span<const std::size_t> DirectedAcyclicGraph::children(std::size_t vertex) const {
  if (vertex >= nvertices()) {
    throw std::out_of_range("Vertex index out of range");
  }
  std::size_t b = csr_.row_ptr().at(vertex);
  std::size_t e = csr_.row_ptr().at(vertex + 1);
  return std::span<const std::size_t>(csr_.col_ind().data() + b, e - b);
}

std::span<const std::size_t> DirectedAcyclicGraph::parents(std::size_t vertex) const {
  if (vertex >= nvertices()) {
    throw std::out_of_range("Vertex index out of range");
  }
  std::size_t b = transpose_.row_ptr().at(vertex);
  std::size_t e = transpose_.row_ptr().at(vertex + 1);
  return std::span<const std::size_t>(transpose_.col_ind().data() + b, e - b);
}

std::size_t DirectedAcyclicGraph::nvertices() const noexcept { return csr_.nvertices(); }
std::size_t DirectedAcyclicGraph::nedges() const noexcept { return csr_.col_ind().size(); }

} // namespace mdgm
