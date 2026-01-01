#include <gtest/gtest.h>

#include <cstddef>
#include <mdgm/graph_storage.hpp>
#include <mdgm/rng.hpp>
#include <mdgm/undirected_graph.hpp>
#include <span>
#include <vector>

namespace mdgm {

TEST(UndirectedGraph, GenerateRegularGraph) {
  std::vector<std::size_t> dims = {4, 4};
  int order = 1;
  UndirectedGraph graph = GenerateRegularGraph(dims, order);

  std:size_t nedges = 0;
  for (std::size_t v = 0; v < graph.nvertices(); ++v) {
    nedges += graph.neighbors(v).size();
  }
  EXPECT_EQ(nedges / 2, graph.nedges());
  EXPECT_EQ(graph.nvertices(), 16);
  EXPECT_EQ(graph.nedges(), 24);
  
  EXPECT_EQ(graph.neighbors(0).size(), 2);
  EXPECT_EQ(graph.neighbors(2).size(), 3);
  EXPECT_EQ(graph.neighbors(5).size(), 4);
  EXPECT_EQ(graph.neighbors(14).size(), 3);
  EXPECT_EQ(graph.neighbors(15).size(), 2);

  std::vector<std::size_t> expected_nbrs_0 = {1, 4};
  std::vector<std::size_t> expected_nbrs_2 = {1, 3, 6};
  std::vector<std::size_t> expected_nbrs_5 = {1, 4, 6, 9};
  std::vector<std::size_t> expected_nbrs_14 = {10, 13, 15};
  std::vector<std::size_t> expected_nbrs_15 = {11, 14};

  auto nbrs_0 = graph.neighbors(0);
  auto nbrs_2 = graph.neighbors(2);
  auto nbrs_5 = graph.neighbors(5);
  auto nbrs_14 = graph.neighbors(14);
  auto nbrs_15 = graph.neighbors(15);

  EXPECT_EQ(std::vector<std::size_t>(nbrs_0.begin(), nbrs_0.end()), expected_nbrs_0);
  EXPECT_EQ(std::vector<std::size_t>(nbrs_2.begin(), nbrs_2.end()), expected_nbrs_2);
  EXPECT_EQ(std::vector<std::size_t>(nbrs_5.begin(), nbrs_5.end()), expected_nbrs_5);
  EXPECT_EQ(std::vector<std::size_t>(nbrs_14.begin(), nbrs_14.end()), expected_nbrs_14);
  EXPECT_EQ(std::vector<std::size_t>(nbrs_15.begin(), nbrs_15.end()), expected_nbrs_15);

  order = 2;
  graph = GenerateRegularGraph(dims, order);

  EXPECT_EQ(graph.nvertices(), 16);
  EXPECT_EQ(graph.nedges(), 42);

  expected_nbrs_0 = {1, 4, 5};
  expected_nbrs_2 = {1, 3, 5, 6, 7};
  expected_nbrs_5 = {0, 1, 2, 4, 6, 8, 9, 10};
  expected_nbrs_14 = {9, 10, 11, 13, 15};
  expected_nbrs_15 = {10, 11, 14};

  nbrs_0 = graph.neighbors(0);
  nbrs_2 = graph.neighbors(2);
  nbrs_5 = graph.neighbors(5);
  nbrs_14 = graph.neighbors(14);
  nbrs_15 = graph.neighbors(15);

  EXPECT_EQ(std::vector<std::size_t>(nbrs_0.begin(), nbrs_0.end()), expected_nbrs_0);
  EXPECT_EQ(std::vector<std::size_t>(nbrs_2.begin(), nbrs_2.end()), expected_nbrs_2);
  EXPECT_EQ(std::vector<std::size_t>(nbrs_5.begin(), nbrs_5.end()), expected_nbrs_5);
  EXPECT_EQ(std::vector<std::size_t>(nbrs_14.begin(), nbrs_14.end()), expected_nbrs_14);
  EXPECT_EQ(std::vector<std::size_t>(nbrs_15.begin(), nbrs_15.end()), expected_nbrs_15);
}

TEST(UndirectedGraph, ValidateUndirected) {
  // Create a simple undirected graph in COO format
  std::vector<std::size_t> row_ind = {0, 0, 1, 1, 2, 2};
  std::vector<std::size_t> col_ind = {1, 2, 0, 2, 0, 1};
  std::vector<double> weights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  GraphCOO coo(3, row_ind, col_ind, weights);

  EXPECT_NO_THROW({
    UndirectedGraph graph(coo);
  });

  // Create a directed graph (missing reverse edges)
  row_ind = {0, 1, 2};
  col_ind = {1, 2, 0};
  weights = {1.0, 1.0, 1.0};
  GraphCOO directed_coo(3, row_ind, col_ind, weights);

  EXPECT_THROW({
    UndirectedGraph graph(directed_coo);
  }, std::invalid_argument);

}

TEST(UndirectedGraph, ValidateConnected) {
  // Create a connected undirected graph in COO format
  std::vector<std::size_t> row_ind = {0, 0, 1, 1, 2, 2};
  std::vector<std::size_t> col_ind = {1, 2, 0, 2, 0, 1};
  std::vector<double> weights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  GraphCOO connected_coo(3, row_ind, col_ind, weights);

  EXPECT_NO_THROW({
    UndirectedGraph graph(connected_coo);
  });

  // Create a disconnected undirected graph
  row_ind = {0, 1};
  col_ind = {1, 0};
  weights = {1.0, 1.0};
  GraphCOO disconnected_coo(3, row_ind, col_ind, weights);

  testing::internal::CaptureStderr();
  EXPECT_NO_THROW({
    UndirectedGraph graph(disconnected_coo);
  });
  std::string output = testing::internal::GetCapturedStderr();
  EXPECT_NE(output.find("warning: graph is not connected"), std::string::npos);
}

TEST(UndirectedGraph, SingleVertexGraph) {
  std::vector<std::size_t> row_ind = {};
  std::vector<std::size_t> col_ind = {};
  std::vector<double> weights = {};
  GraphCOO coo(1, row_ind, col_ind, weights);

  EXPECT_NO_THROW({
    UndirectedGraph graph(coo);
  });

  UndirectedGraph graph(coo);
  EXPECT_EQ(graph.nvertices(), 1);
  EXPECT_EQ(graph.nedges(), 0);
  EXPECT_EQ(graph.neighbors(0).size(), 0);
}

TEST(UndirectedGraph, CompleteGraph) {
  // Complete graph K4: every vertex connected to every other vertex
  std::vector<std::size_t> row_ind = {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3};
  std::vector<std::size_t> col_ind = {1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2};
  std::vector<double> weights(12, 1.0);
  GraphCOO coo(4, row_ind, col_ind, weights);

  EXPECT_NO_THROW({
    UndirectedGraph graph(coo);
  });

  UndirectedGraph graph(coo);
  EXPECT_EQ(graph.nvertices(), 4);
  EXPECT_EQ(graph.nedges(), 6);
  for (std::size_t v = 0; v < 4; ++v) {
    EXPECT_EQ(graph.neighbors(v).size(), 3);
  }
}

TEST(UndirectedGraph, SpanningTreeSamplers) {
  // Test spanning tree for larger graph (4-cycle)
  std::vector<std::size_t> row_ind = {0, 1, 2, 3, 1, 2, 3, 0};
  std::vector<std::size_t> col_ind = {1, 0, 1, 2, 2, 3, 0, 3};
  std::vector<double> weights(8, 1.0);
  GraphCOO coo(4, row_ind, col_ind, weights);
  UndirectedGraph graph(coo);

  RNG rng(42);
  GraphCOO spanning_tree = graph.SampleSpanningTree(rng, kWilson);

  // Spanning tree must have n-1 edges
  EXPECT_EQ(spanning_tree.nvertices(), 4);
  EXPECT_EQ(spanning_tree.row_ind().size(), 3);
  EXPECT_EQ(spanning_tree.col_ind().size(), 3);

  spanning_tree = graph.SampleSpanningTree(rng, kAldousBroder);

  // Spanning tree must have n-1 edges
  EXPECT_EQ(spanning_tree.nvertices(), 4);
  EXPECT_EQ(spanning_tree.row_ind().size(), 3);
  EXPECT_EQ(spanning_tree.col_ind().size(), 3);
}

TEST(UndirectedGraph, DeterministicRandomness) {
  UndirectedGraph graph = GenerateRegularGraph({3, 3}, 1);

  RNG rng1(42);
  GraphCOO tree1 = graph.SampleSpanningTree(rng1, kWilson);

  RNG rng2(42);
  GraphCOO tree2 = graph.SampleSpanningTree(rng2, kWilson);

  // Same seed should produce same spanning tree
  EXPECT_EQ(tree1.row_ind(), tree2.row_ind());
  EXPECT_EQ(tree1.col_ind(), tree2.col_ind());
}

}  // namespace mdgm
