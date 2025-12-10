#include <algorithm>
#include <cstddef>
#include <mdgm/graph_storage.hpp>
#include <span>
#include <vector>

namespace mdgm {

GraphCSR::GraphCSR(const std::size_t nvertices, 
                   const std::vector<std::size_t>& row_ptr,
                   const std::vector<std::size_t>& col_ind,
                   const std::vector<double>& weights)
    : nvertices_(nvertices),
      row_ptr_(row_ptr),
      col_ind_(col_ind),
      weights_(weights) {
        if (row_ptr_.size() != nvertices_ + 1) {
          throw std::invalid_argument("row_ptr size must be nvertices + 1");
        }
        if (col_ind_.size() != weights_.size()) {
          throw std::invalid_argument("col_ind and weights must have the same length");
        }
        if (row_ptr_.back() != col_ind_.size()) {
          throw std::invalid_argument("row_ptr last element must equal number of edges");
        } 
      }

GraphCSR::GraphCSR(const GraphCOO& coo)
    : nvertices_(coo.nvertices()),
      row_ptr_(nvertices_ + 1, 0),
      col_ind_(coo.col_ind().size(), 0),
      weights_(coo.weights().size(), 0) {
  BuildFromCOO_(coo);
  SortAndDeduplicate_();
}

void GraphCSR::BuildFromCOO_(const GraphCOO& coo) {
  for (const auto& r : coo.row_ind()) {
    ++row_ptr_[r + 1];
  }

  for (std::size_t i = 1; i <= nvertices_; ++i) {
    row_ptr_[i] += row_ptr_[i - 1];
  }

  std::vector<std::size_t> counter = row_ptr_;

  for (std::size_t i = 0; i < coo.col_ind().size(); ++i) {
    std::size_t row = coo.row_ind().at(i);
    std::size_t dest = counter[row]++;
    col_ind_[dest] = coo.col_ind().at(i);
    weights_[dest] = coo.weights().at(i);
  }
}

void GraphCSR::SortAndDeduplicate_() {
  // NOTE: Duplicate edges are discarded; only the first occruence (lowest column after sort) is
  // kept. Weights of discarded edges are ignored. 
  
  // Store original row_ptr_ for correct reference of col_ind_ and weights_ during processing
  std::vector<std::size_t> orig_row_ptr_ = row_ptr_;
  std::vector<std::pair<std::size_t, double>> tmp;
  for (std::size_t i = 0; i < nvertices_; ++i) {
    std::size_t b = orig_row_ptr_[i];
    std::size_t e = orig_row_ptr_[i + 1];

    if (e - b > 1) {
      tmp.clear();
      tmp.reserve(e - b);
      for (std::size_t j = b; j < e; ++j) {
        tmp.emplace_back(col_ind_[j], weights_[j]);
      }

      // sorts by column index
      std::sort(tmp.begin(), tmp.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });

      std::size_t k = row_ptr_[i];  // position to write deduplicated entries
      for (std::size_t j = 0; j < tmp.size(); ++j) {
        if (j == 0 || tmp[j].first != tmp[j - 1].first) {
          col_ind_[k] = tmp[j].first;
          weights_[k] = tmp[j].second;
          ++k;
        }
      }
      row_ptr_[i + 1] = k;  // Update row_ptr_ to new end position
    } else {
      row_ptr_[i + 1] = row_ptr_[i] + (e - b);  // adjust for previous duplicates
    }
  }

  std::size_t new_size = row_ptr_[nvertices_];
  col_ind_.resize(new_size);
  weights_.resize(new_size);
}

std::span<const std::size_t> GraphCSR::adjacent(std::size_t vertex) const {
  if (vertex >= nvertices_) {
    return {};
  }
  std::size_t b = row_ptr_[vertex];
  std::size_t e = row_ptr_[vertex + 1];
  return std::span<const std::size_t>(col_ind_.data() + b, e - b);
}

std::span<const double> GraphCSR::adjacent_weights(std::size_t vertex) const {
  if (vertex >= nvertices_) {
    return {};
  }
  std::size_t b = row_ptr_[vertex];
  std::size_t e = row_ptr_[vertex + 1];
  return std::span<const double>(weights_.data() + b, e - b);
}

}  // namespace mdgm
