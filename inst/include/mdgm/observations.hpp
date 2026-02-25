#pragma once

#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace mdgm {

// Compressed storage for per-vertex observations.
// Same pattern as CSR row_ptr: obs_ptr[i]..obs_ptr[i+1] indexes into obs_data.
struct Observations {
  std::vector<int> data;
  std::vector<std::size_t> ptr;  // size = nvertices + 1

  std::size_t nvertices() const { return ptr.size() > 0 ? ptr.size() - 1 : 0; }

  std::span<const int> operator[](std::size_t vertex) const {
    return {data.data() + ptr[vertex], ptr[vertex + 1] - ptr[vertex]};
  }

  bool empty(std::size_t vertex) const {
    return ptr[vertex] == ptr[vertex + 1];
  }
};

}  // namespace mdgm
