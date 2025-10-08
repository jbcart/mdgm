#include <cstddef>
#include <vector>

class Graph {
 public:
  Graph(std::size_t num_rows, std::size_t num_cols, int order);
  ~Graph();

 private:
  std::vector<std::size_t> row_ptr {};
  std::vector<std::size_t> col_ind {};
  std::vector<double> values {};
};