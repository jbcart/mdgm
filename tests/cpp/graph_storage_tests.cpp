#include <gtest/gtest.h>

#include <cstddef>
#include <mdgm/graph_storage.hpp>
#include <span>
#include <vector>

namespace mdgm {

TEST(GraphStorage, COOConstruction) {
  std::size_t nvertices = 3;
  std::vector<std::size_t> row_ind = {0, 1, 2};
  std::vector<std::size_t> col_ind = {1, 2};
  std::vector<double> weights_eq = {1.0, 1.0};
  std::vector<double> weights = {1.0, 2.0, 3.0};

  EXPECT_THROW({ GraphCOO coo_graph(nvertices, row_ind, col_ind); }, std::invalid_argument);
  col_ind.push_back(0);

  EXPECT_THROW(
      { GraphCOO coo_graph(nvertices, row_ind, col_ind, weights_eq); }, std::invalid_argument);
  weights_eq.push_back(1.0);

  GraphCOO coo_graph(nvertices, row_ind, col_ind);

  EXPECT_EQ(coo_graph.nvertices(), nvertices);
  EXPECT_EQ(coo_graph.row_ind(), row_ind);
  EXPECT_EQ(coo_graph.col_ind(), col_ind);
  EXPECT_EQ(coo_graph.weights(), weights_eq);

  GraphCOO coo_graph_w(nvertices, row_ind, col_ind, weights);
  EXPECT_EQ(coo_graph_w.weights(), weights);
}

TEST(GraphStorage, COOtoCSR) {
  // Create a COO graph with some duplicate edges
  std::size_t nvertices = 4;
  std::vector<std::size_t> row_ind = {1, 2, 1, 1, 0, 0, 2};
  std::vector<std::size_t> col_ind = {3, 3, 2, 2, 3, 1, 3};
  std::vector<double> weights = {4.0, 5.0, 3.0, 3.0, 2.0, 1.0, 5.0};

  GraphCOO coo_graph(nvertices, row_ind, col_ind, weights);
  GraphCSR csr_graph(coo_graph);

  // Check CSR structure - sorted and deduplicated
  std::vector<std::size_t> expected_row_ptr = {0, 2, 4, 5, 5};
  std::vector<std::size_t> expected_col_ind = {1, 3, 2, 3, 3};
  std::vector<double> expected_weights = {1.0, 2.0, 3.0, 4.0, 5.0};

  EXPECT_EQ(csr_graph.row_ptr(), expected_row_ptr);
  EXPECT_EQ(csr_graph.col_ind(), expected_col_ind);
  EXPECT_EQ(csr_graph.weights(), expected_weights);

  // Check adjacency retrieval
  auto adj_vertex_0 = csr_graph.adjacent(0);
  std::vector<std::size_t> expected_adj_0 = {1, 3};
  EXPECT_EQ(std::vector<std::size_t>(adj_vertex_0.begin(), adj_vertex_0.end()), expected_adj_0);

  auto adj_weights_0 = csr_graph.adjacent_weights(0);
  std::vector<double> expected_adj_weights_0 = {1.0, 2.0};
  EXPECT_EQ(std::vector<double>(adj_weights_0.begin(), adj_weights_0.end()),
            expected_adj_weights_0);
}

TEST(GraphStorage, CSRConstruction) {
  std::size_t nvertices = 3;
  std::vector<std::size_t> row_ptr = {0, 2, 3};
  std::vector<std::size_t> col_ind = {1, 2, 0};
  std::vector<double> weights = {1.0, 2.0, 3.0};
  std::vector<double> weights_bad = {1.0, 2.0};

  EXPECT_THROW(
      { GraphCSR csr_graph(nvertices, row_ptr, col_ind, weights); }, std::invalid_argument);
  row_ptr.push_back(4);
  EXPECT_THROW(
      { GraphCSR csr_graph(nvertices, row_ptr, col_ind, weights); }, std::invalid_argument);
  row_ptr.pop_back();
  row_ptr.push_back(3);  // correct size now
  EXPECT_THROW(
      { GraphCSR csr_graph(nvertices, row_ptr, col_ind, weights_bad); }, std::invalid_argument);

  GraphCSR csr_graph(nvertices, row_ptr, col_ind, weights);

  EXPECT_EQ(csr_graph.nvertices(), nvertices);
  EXPECT_EQ(csr_graph.row_ptr(), row_ptr);
  EXPECT_EQ(csr_graph.col_ind(), col_ind);
  EXPECT_EQ(csr_graph.weights(), weights);
  EXPECT_EQ(csr_graph.row_ptr().back(), col_ind.size());
}

}  // namespace mdgm
