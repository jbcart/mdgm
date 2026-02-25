#include <gtest/gtest.h>

#include <cstddef>
#include <mdgm/directed_acyclic_graph.hpp>
#include <mdgm/graph_storage.hpp>
#include <span>
#include <vector>

namespace mdgm {

// Helper: convert span to vector for EXPECT_EQ comparisons
static std::vector<std::size_t> ToVec(std::span<const std::size_t> s) {
  return {s.begin(), s.end()};
}

TEST(DirectedAcyclicGraph, BasicConstruction) {
  // Linear chain: 0 → 1 → 2
  std::vector<std::size_t> row_ind = {0, 1};
  std::vector<std::size_t> col_ind = {1, 2};
  GraphCOO coo(3, row_ind, col_ind);

  EXPECT_NO_THROW({ DirectedAcyclicGraph dag(coo, DagType::kSpanningTree); });

  DirectedAcyclicGraph dag(coo, DagType::kSpanningTree);
  EXPECT_EQ(dag.nvertices(), 3);
  EXPECT_EQ(dag.nedges(), 2);
}

TEST(DirectedAcyclicGraph, ChildrenQuery) {
  // Linear chain: 0 → 1 → 2
  std::vector<std::size_t> row_ind = {0, 1};
  std::vector<std::size_t> col_ind = {1, 2};
  GraphCOO coo(3, row_ind, col_ind);
  DirectedAcyclicGraph dag(coo, DagType::kSpanningTree);

  EXPECT_EQ(ToVec(dag.children(0)), (std::vector<std::size_t>{1}));
  EXPECT_EQ(ToVec(dag.children(1)), (std::vector<std::size_t>{2}));
  EXPECT_EQ(ToVec(dag.children(2)), (std::vector<std::size_t>{}));
}

TEST(DirectedAcyclicGraph, ParentsQuery) {
  // Linear chain: 0 → 1 → 2
  std::vector<std::size_t> row_ind = {0, 1};
  std::vector<std::size_t> col_ind = {1, 2};
  GraphCOO coo(3, row_ind, col_ind);
  DirectedAcyclicGraph dag(coo, DagType::kSpanningTree);

  EXPECT_EQ(ToVec(dag.parents(0)), (std::vector<std::size_t>{}));
  EXPECT_EQ(ToVec(dag.parents(1)), (std::vector<std::size_t>{0}));
  EXPECT_EQ(ToVec(dag.parents(2)), (std::vector<std::size_t>{1}));
}

TEST(DirectedAcyclicGraph, ChildrenParentsConsistency) {
  // Diamond DAG: 0 → {1, 2}, 1 → 3, 2 → 3
  std::vector<std::size_t> row_ind = {0, 0, 1, 2};
  std::vector<std::size_t> col_ind = {1, 2, 3, 3};
  GraphCOO coo(4, row_ind, col_ind);
  DirectedAcyclicGraph dag(coo, DagType::kAcyclicOrientation);

  // For every u, for every v in children(u): u must be in parents(v)
  for (std::size_t u = 0; u < dag.nvertices(); ++u) {
    for (std::size_t v : dag.children(u)) {
      auto p = ToVec(dag.parents(v));
      EXPECT_NE(std::find(p.begin(), p.end(), u), p.end())
          << "u=" << u << " is in children but not in parents of v=" << v;
    }
  }

  // For every v, for every u in parents(v): v must be in children(u)
  for (std::size_t v = 0; v < dag.nvertices(); ++v) {
    for (std::size_t u : dag.parents(v)) {
      auto c = ToVec(dag.children(u));
      EXPECT_NE(std::find(c.begin(), c.end(), v), c.end())
          << "u=" << u << " is in parents but v=" << v << " not in children";
    }
  }
}

TEST(DirectedAcyclicGraph, OutOfRangeAccess) {
  std::vector<std::size_t> row_ind = {0};
  std::vector<std::size_t> col_ind = {1};
  GraphCOO coo(2, row_ind, col_ind);
  DirectedAcyclicGraph dag(coo, DagType::kSpanningTree);

  EXPECT_THROW(dag.children(2), std::out_of_range);
  EXPECT_THROW(dag.parents(2), std::out_of_range);
  EXPECT_THROW(dag.children(dag.nvertices()), std::out_of_range);
  EXPECT_THROW(dag.parents(dag.nvertices()), std::out_of_range);
}

TEST(DirectedAcyclicGraph, SingleVertexGraph) {
  GraphCOO coo(1, {}, {});

  EXPECT_NO_THROW({ DirectedAcyclicGraph dag(coo, DagType::kRooted); });

  DirectedAcyclicGraph dag(coo, DagType::kRooted);
  EXPECT_EQ(dag.nvertices(), 1);
  EXPECT_EQ(dag.nedges(), 0);
  EXPECT_EQ(ToVec(dag.children(0)), (std::vector<std::size_t>{}));
  EXPECT_EQ(ToVec(dag.parents(0)), (std::vector<std::size_t>{}));
}

TEST(DirectedAcyclicGraph, BranchingDAG) {
  // Diamond DAG: 0 → {1, 2}, 1 → 3, 2 → 3
  std::vector<std::size_t> row_ind = {0, 0, 1, 2};
  std::vector<std::size_t> col_ind = {1, 2, 3, 3};
  GraphCOO coo(4, row_ind, col_ind);
  DirectedAcyclicGraph dag(coo, DagType::kAcyclicOrientation);

  EXPECT_EQ(dag.nvertices(), 4);
  EXPECT_EQ(dag.nedges(), 4);

  EXPECT_EQ(ToVec(dag.children(0)), (std::vector<std::size_t>{1, 2}));
  EXPECT_EQ(ToVec(dag.children(1)), (std::vector<std::size_t>{3}));
  EXPECT_EQ(ToVec(dag.children(2)), (std::vector<std::size_t>{3}));
  EXPECT_EQ(ToVec(dag.children(3)), (std::vector<std::size_t>{}));

  EXPECT_EQ(ToVec(dag.parents(0)), (std::vector<std::size_t>{}));
  EXPECT_EQ(ToVec(dag.parents(1)), (std::vector<std::size_t>{0}));
  EXPECT_EQ(ToVec(dag.parents(2)), (std::vector<std::size_t>{0}));
  EXPECT_EQ(ToVec(dag.parents(3)), (std::vector<std::size_t>{1, 2}));
}

TEST(DirectedAcyclicGraph, CopySemantics) {
  std::vector<std::size_t> row_ind = {0, 0, 1};
  std::vector<std::size_t> col_ind = {1, 2, 2};
  GraphCOO coo(3, row_ind, col_ind);
  DirectedAcyclicGraph dag(coo, DagType::kSpanningTree);

  DirectedAcyclicGraph dag_copy = dag;  // copy constructor
  EXPECT_EQ(dag_copy.nvertices(), dag.nvertices());
  EXPECT_EQ(dag_copy.nedges(), dag.nedges());
  for (std::size_t v = 0; v < dag.nvertices(); ++v) {
    EXPECT_EQ(ToVec(dag_copy.children(v)), ToVec(dag.children(v)));
    EXPECT_EQ(ToVec(dag_copy.parents(v)), ToVec(dag.parents(v)));
  }
}

}  // namespace mdgm
