#pragma once

#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace mdgm {

class GraphCOO {
 public:
  GraphCOO(std::size_t nvertices, const std::vector<std::size_t>& row_ind,
           const std::vector<std::size_t>& col_ind,
           const std::vector<double>& weights)
      : nvertices_(nvertices),
        row_ind_(row_ind),
        col_ind_(col_ind),
        weights_(weights) {
    if (row_ind.size() != col_ind.size() || row_ind.size() != weights.size()) {
      throw std::invalid_argument("COO arrays must have the same length");
    }
  }

  GraphCOO(std::size_t nvertices, const std::vector<std::size_t>& row_ind,
           const std::vector<std::size_t>& col_ind)
      : GraphCOO(nvertices, row_ind, col_ind,
                 std::vector<double>(row_ind.size(), 1.0)) {}

  ~GraphCOO() = default;
  GraphCOO(const GraphCOO&) = default;
  GraphCOO& operator=(const GraphCOO&) = default;

  void add_edge(std::size_t row, std::size_t col, double weight) {
    row_ind_.push_back(row);
    col_ind_.push_back(col);
    weights_.push_back(weight);
  }

  std::size_t nvertices() const noexcept { return nvertices_; }
  const std::vector<std::size_t>& row_ind() const noexcept { return row_ind_; }
  const std::vector<std::size_t>& col_ind() const noexcept { return col_ind_; }
  const std::vector<double>& weights() const noexcept { return weights_; }

 private:
  std::size_t nvertices_;
  std::vector<std::size_t> row_ind_;
  std::vector<std::size_t> col_ind_;
  std::vector<double> weights_;
};

class GraphCSR {
 public:
  explicit GraphCSR(const GraphCOO& coo);

  ~GraphCSR() = default;
  GraphCSR(const GraphCSR&) = default;
  GraphCSR& operator=(const GraphCSR&) = default;

  // CSR accessors
  const std::vector<std::size_t>& row_ptr() const noexcept { return row_ptr_; }
  const std::vector<std::size_t>& col_ind() const noexcept { return col_ind_; }
  const std::vector<double>& weights() const noexcept { return weights_; }
  const std::size_t nvertices() const noexcept { return nvertices_; }

  // Get adjacent vertices and weights for a given vertex
  std::span<const std::size_t> adjacent(std::size_t vertex) const;
  std::span<const double> adjacent_weights(std::size_t vertex) const;

 private:
  void BuildFromCOO_(const GraphCOO& coo);
  void SortAndDeduplicate_();

  std::size_t nvertices_;
  std::vector<std::size_t> row_ptr_;
  std::vector<std::size_t> col_ind_;
  std::vector<double> weights_;
};

}  // namespace mdgm
