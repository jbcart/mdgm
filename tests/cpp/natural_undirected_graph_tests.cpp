#include <gtest/gtest.h>

#include <cstddef>
#include <mdgm/graph_storage.hpp>
#include <mdgm/rng.hpp>
#include <mdgm/natural_undirected_graph.hpp>
#include <mdgm/directed_acyclic_graph.hpp>
#include <span>
#include <vector>

namespace mdgm {

TEST(NaturalUndirectedGraph, GenerateRegularGraph) {
  std::vector<std::size_t> dims = {4, 4};
  int order = 1;
  NaturalUndirectedGraph graph = GenerateRegularGraph(dims, order);

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

TEST(NaturalUndirectedGraph, ValidateUndirected) {
  // Create a simple undirected graph in COO format
  std::vector<std::size_t> row_ind = {0, 0, 1, 1, 2, 2};
  std::vector<std::size_t> col_ind = {1, 2, 0, 2, 0, 1};
  std::vector<double> weights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  GraphCOO coo(3, row_ind, col_ind, weights);

  EXPECT_NO_THROW({
    NaturalUndirectedGraph graph(coo);
  });

  // Create a directed graph (missing reverse edges)
  row_ind = {0, 1, 2};
  col_ind = {1, 2, 0};
  weights = {1.0, 1.0, 1.0};
  GraphCOO directed_coo(3, row_ind, col_ind, weights);

  EXPECT_THROW({
    NaturalUndirectedGraph graph(directed_coo);
  }, std::invalid_argument);

}

TEST(NaturalUndirectedGraph, ValidateConnected) {
  // Create a connected undirected graph in COO format
  std::vector<std::size_t> row_ind = {0, 0, 1, 1, 2, 2};
  std::vector<std::size_t> col_ind = {1, 2, 0, 2, 0, 1};
  std::vector<double> weights = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  GraphCOO connected_coo(3, row_ind, col_ind, weights);

  EXPECT_NO_THROW({
    NaturalUndirectedGraph graph(connected_coo);
  });

  // Create a disconnected undirected graph
  row_ind = {0, 1};
  col_ind = {1, 0};
  weights = {1.0, 1.0};
  GraphCOO disconnected_coo(3, row_ind, col_ind, weights);

  // Construction succeeds but spanning tree sampling should throw
  NaturalUndirectedGraph graph(disconnected_coo);
  EXPECT_EQ(graph.nvertices(), 3);
  RNG rng(42);
  EXPECT_THROW(graph.SampleSpanningTree(rng, SpanningTreeMethod::kWilson), std::runtime_error);
}

TEST(NaturalUndirectedGraph, SingleVertexGraph) {
  std::vector<std::size_t> row_ind = {};
  std::vector<std::size_t> col_ind = {};
  std::vector<double> weights = {};
  GraphCOO coo(1, row_ind, col_ind, weights);

  EXPECT_NO_THROW({
    NaturalUndirectedGraph graph(coo);
  });

  NaturalUndirectedGraph graph(coo);
  EXPECT_EQ(graph.nvertices(), 1);
  EXPECT_EQ(graph.nedges(), 0);
  EXPECT_EQ(graph.neighbors(0).size(), 0);
}

TEST(NaturalUndirectedGraph, CompleteGraph) {
  // Complete graph K4: every vertex connected to every other vertex
  std::vector<std::size_t> row_ind = {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3};
  std::vector<std::size_t> col_ind = {1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2};
  std::vector<double> weights(12, 1.0);
  GraphCOO coo(4, row_ind, col_ind, weights);

  EXPECT_NO_THROW({
    NaturalUndirectedGraph graph(coo);
  });

  NaturalUndirectedGraph graph(coo);
  EXPECT_EQ(graph.nvertices(), 4);
  EXPECT_EQ(graph.nedges(), 6);
  for (std::size_t v = 0; v < 4; ++v) {
    EXPECT_EQ(graph.neighbors(v).size(), 3);
  }
}

TEST(NaturalUndirectedGraph, SpanningTreeSamplers) {
  // Test spanning tree for larger graph (4-cycle)
  std::vector<std::size_t> row_ind = {0, 1, 2, 3, 1, 2, 3, 0};
  std::vector<std::size_t> col_ind = {1, 0, 1, 2, 2, 3, 0, 3};
  std::vector<double> weights(8, 1.0);
  GraphCOO coo(4, row_ind, col_ind, weights);
  NaturalUndirectedGraph graph(coo);

  RNG rng(42);
  DirectedAcyclicGraph spanning_tree = graph.SampleSpanningTree(rng, SpanningTreeMethod::kWilson);

  // Spanning tree must have n-1 edges
  EXPECT_EQ(spanning_tree.nvertices(), 4);
  EXPECT_EQ(spanning_tree.nedges(), 3);
  for (std::size_t v = 0; v < spanning_tree.nvertices(); ++v) {
    EXPECT_LE(spanning_tree.parents(v).size(), 1);
  }

  spanning_tree = graph.SampleSpanningTree(rng, SpanningTreeMethod::kAldousBroder);

  // Spanning tree must have n-1 edges
  EXPECT_EQ(spanning_tree.nvertices(), 4);
  EXPECT_EQ(spanning_tree.nedges(), 3);
  for (std::size_t v = 0; v < spanning_tree.nvertices(); ++v) {
    EXPECT_LE(spanning_tree.parents(v).size(), 1);
  }
}

TEST(NaturalUndirectedGraph, DeterministicRandomness) {
  NaturalUndirectedGraph graph = GenerateRegularGraph({3, 3}, 1);

  RNG rng1(42);
  DirectedAcyclicGraph tree1 = graph.SampleSpanningTree(rng1, SpanningTreeMethod::kWilson);
  DirectedAcyclicGraph ao1 = graph.SampleAcyclicOrientation(rng1);

  RNG rng2(42);
  DirectedAcyclicGraph tree2 = graph.SampleSpanningTree(rng2, SpanningTreeMethod::kWilson);
  DirectedAcyclicGraph ao2 = graph.SampleAcyclicOrientation(rng2);

  // Same seed should produce same spanning tree
  for (std::size_t i = 0; i < tree1.nvertices(); ++i) {
    for (std::size_t j = 0; j < tree1.children(i).size(); ++j) {
      EXPECT_EQ(tree1.children(i)[j], tree2.children(i)[j]);
    }
    for (std::size_t j = 0; j < tree1.parents(i).size(); ++j) {
      EXPECT_EQ(tree1.parents(i)[j], tree2.parents(i)[j]);
    }
    for (std::size_t j = 0; j < ao1.children(i).size(); ++j) {
      EXPECT_EQ(ao1.children(i)[j], ao2.children(i)[j]);
    }
    for (std::size_t j = 0; j < ao1.parents(i).size(); ++j) {
      EXPECT_EQ(ao1.parents(i)[j], ao2.parents(i)[j]);
    }
  }
} 

TEST(NaturalUndirectedGraph, AcyclicOrientation) {
  NaturalUndirectedGraph graph = GenerateRegularGraph({3, 3}, 1);
  RNG rng(42);
  DirectedAcyclicGraph acyclic_orientation = graph.SampleAcyclicOrientation(rng);

  // Check that the acyclic orientation has the same number of vertices and edges
  EXPECT_EQ(acyclic_orientation.nvertices(), graph.nvertices());
  EXPECT_EQ(acyclic_orientation.nedges(), graph.nedges());
}

}  // namespace mdgm
