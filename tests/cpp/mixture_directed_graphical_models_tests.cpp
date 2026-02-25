#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <mdgm/directed_acyclic_graph.hpp>
#include <mdgm/graph_storage.hpp>
#include <mdgm/mixture_directed_graphical_models.hpp>
#include <mdgm/natural_undirected_graph.hpp>
#include <mdgm/rng.hpp>
#include <vector>

namespace mdgm {

// Helper to create a small triangle graph (K3): 0-1-2-0
static NaturalUndirectedGraph MakeTriangle() {
  std::vector<std::size_t> row = {0, 0, 1, 1, 2, 2};
  std::vector<std::size_t> col = {1, 2, 0, 2, 0, 1};
  return NaturalUndirectedGraph(GraphCOO(3, row, col));
}

// Helper to create a 3x3 grid graph (order 1)
static NaturalUndirectedGraph MakeGrid3x3() {
  return GenerateRegularGraph({3, 3}, 1);
}

TEST(MixtureDirectedGraphicalModels, ConstructionSpanningTree) {
  auto nug = MakeTriangle();
  EXPECT_NO_THROW({
    MixtureDirectedGraphicalModels mdgm(nug, DagType::kSpanningTree, 2);
  });

  MixtureDirectedGraphicalModels mdgm(nug, DagType::kSpanningTree, 2);
  EXPECT_EQ(mdgm.nvertices(), 3);
  EXPECT_EQ(mdgm.ncolors(), 2);
}

TEST(MixtureDirectedGraphicalModels, ConstructionAcyclicOrientation) {
  auto nug = MakeTriangle();
  EXPECT_NO_THROW({
    MixtureDirectedGraphicalModels mdgm(nug, DagType::kAcyclicOrientation, 2);
  });
}

TEST(MixtureDirectedGraphicalModels, ZFullConditionalSpanningTree) {
  // Create a chain graph: 0-1-2
  std::vector<std::size_t> row = {0, 1, 1, 2};
  std::vector<std::size_t> col = {1, 0, 2, 1};
  NaturalUndirectedGraph nug(GraphCOO(3, row, col));

  MixtureDirectedGraphicalModels mdgm(nug, DagType::kSpanningTree, 2);

  // For spanning tree on a chain, the DAG IS the chain (only spanning tree)
  // With z = {0, 0, 1} and psi = 1.0:
  std::vector<int> z = {0, 0, 1};
  double psi = 1.0;

  auto probs = mdgm.ZFullConditional(z, 1, psi);

  // Vertex 1 has neighbors 0 and 2 in the tree
  // z[0] = 0, z[2] = 1
  // p(z_1 = 0) proportional to exp(1.0 * 1 + 0) = e  (matches z[0]=0)
  // p(z_1 = 1) proportional to exp(1.0 * 1 + 0) = e  (matches z[2]=1)
  // So probabilities should be 0.5, 0.5
  EXPECT_EQ(probs.size(), 2);
  EXPECT_NEAR(probs[0], 0.5, 1e-10);
  EXPECT_NEAR(probs[1], 0.5, 1e-10);

  // With z = {0, 0, 0} and psi = 1.0:
  z = {0, 0, 0};
  probs = mdgm.ZFullConditional(z, 1, psi);
  // p(z_1 = 0) proportional to exp(1.0 * 2) = e^2  (matches both neighbors)
  // p(z_1 = 1) proportional to exp(1.0 * 0) = 1     (matches neither)
  // p(0) = e^2 / (e^2 + 1)
  double e2 = std::exp(2.0);
  EXPECT_NEAR(probs[0], e2 / (e2 + 1.0), 1e-10);
  EXPECT_NEAR(probs[1], 1.0 / (e2 + 1.0), 1e-10);
}

TEST(MixtureDirectedGraphicalModels, ZFullConditionalSumsToOne) {
  auto nug = MakeGrid3x3();
  MixtureDirectedGraphicalModels mdgm(nug, DagType::kSpanningTree, 2);

  std::vector<int> z = {0, 1, 0, 1, 0, 1, 0, 1, 0};
  double psi = 0.5;

  for (std::size_t v = 0; v < 9; ++v) {
    auto probs = mdgm.ZFullConditional(z, v, psi);
    double sum = 0.0;
    for (double p : probs) sum += p;
    EXPECT_NEAR(sum, 1.0, 1e-10) << "vertex=" << v;
    for (double p : probs) {
      EXPECT_GE(p, 0.0) << "vertex=" << v;
    }
  }
}

TEST(MixtureDirectedGraphicalModels, ZFullConditionalAOSumsToOne) {
  auto nug = MakeGrid3x3();
  MixtureDirectedGraphicalModels mdgm(nug, DagType::kAcyclicOrientation, 2);

  std::vector<int> z = {0, 1, 0, 1, 0, 1, 0, 1, 0};
  double psi = 0.5;

  for (std::size_t v = 0; v < 9; ++v) {
    auto probs = mdgm.ZFullConditional(z, v, psi);
    double sum = 0.0;
    for (double p : probs) sum += p;
    EXPECT_NEAR(sum, 1.0, 1e-10) << "vertex=" << v;
    for (double p : probs) {
      EXPECT_GE(p, 0.0) << "vertex=" << v;
    }
  }
}

TEST(MixtureDirectedGraphicalModels, LogLikelihoodFinite) {
  auto nug = MakeGrid3x3();
  MixtureDirectedGraphicalModels mdgm(nug, DagType::kSpanningTree, 2);

  std::vector<int> z = {0, 1, 0, 1, 0, 1, 0, 1, 0};
  double ll = mdgm.LogLikelihood(z, 0.5);
  EXPECT_TRUE(std::isfinite(ll));
  EXPECT_LT(ll, 0.0);  // log-likelihood should be negative
}

TEST(MixtureDirectedGraphicalModels, LogLikelihoodZeroPsi) {
  // When psi = 0, all conditionals are uniform: p(z_i|...) = 1/n_colors
  auto nug = MakeGrid3x3();
  MixtureDirectedGraphicalModels mdgm(nug, DagType::kSpanningTree, 2);

  std::vector<int> z = {0, 1, 0, 1, 0, 1, 0, 1, 0};
  double ll = mdgm.LogLikelihood(z, 0.0);
  // With psi=0, alpha all 0: each conditional is 1/2
  // ll = 9 * log(1/2) = -9 * log(2)
  EXPECT_NEAR(ll, -9.0 * std::log(2.0), 1e-10);
}

TEST(MixtureDirectedGraphicalModels, UpdateGraphSpanningTree) {
  auto nug = MakeGrid3x3();
  MixtureDirectedGraphicalModels mdgm(nug, DagType::kSpanningTree, 2);
  RNG rng(42);

  std::vector<int> z = {0, 0, 0, 1, 1, 1, 0, 0, 0};
  double psi = 1.0;

  // UpdateGraph should not throw and should produce a valid spanning tree
  EXPECT_NO_THROW(mdgm.UpdateGraph(z, psi, rng));

  // Verify the current DAG is a valid spanning tree (n-1 edges)
  const auto& dag = mdgm.current_dag();
  EXPECT_EQ(dag.nvertices(), 9);
  EXPECT_EQ(dag.nedges(), 8);
}

TEST(MixtureDirectedGraphicalModels, UpdateGraphAO) {
  auto nug = MakeGrid3x3();
  MixtureDirectedGraphicalModels mdgm(nug, DagType::kAcyclicOrientation, 2);
  RNG rng(42);

  std::vector<int> z = {0, 0, 0, 1, 1, 1, 0, 0, 0};
  double psi = 0.5;

  EXPECT_NO_THROW(mdgm.UpdateGraph(z, psi, rng));

  // AO should have same number of edges as NUG
  const auto& dag = mdgm.current_dag();
  EXPECT_EQ(dag.nvertices(), 9);
  EXPECT_EQ(dag.nedges(), nug.nedges());
}

TEST(MixtureDirectedGraphicalModels, StoreSample) {
  auto nug = MakeGrid3x3();
  MixtureDirectedGraphicalModels mdgm(nug, DagType::kSpanningTree, 2);

  std::size_t n = 9;
  std::size_t J = 3;
  std::vector<std::size_t> dag_data(n * J, 0);

  mdgm.StoreSample(dag_data, 0, n);

  // Each vertex should have at most 1 parent (spanning tree)
  std::size_t root_count = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (dag_data[i] == SIZE_MAX) {
      ++root_count;
    } else {
      EXPECT_LT(dag_data[i], n) << "vertex=" << i;
    }
  }
  // Spanning tree has exactly one root
  EXPECT_EQ(root_count, 1);
}

}  // namespace mdgm
