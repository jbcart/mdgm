#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <mdgm/directed_acyclic_graph.hpp>
#include <mdgm/emission.hpp>
#include <mdgm/graph_storage.hpp>
#include <mdgm/mixture_directed_graphical_models.hpp>
#include <mdgm/model.hpp>
#include <mdgm/natural_undirected_graph.hpp>
#include <mdgm/observations.hpp>
#include <mdgm/rng.hpp>
#include <memory>
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

// --- Tests for n_colors > 2 ---

TEST(MixtureDirectedGraphicalModels, ThreeColorZFullConditionalSumsToOne) {
  auto nug = MakeGrid3x3();
  MixtureDirectedGraphicalModels mdgm(nug, DagType::kSpanningTree, 3);

  std::vector<int> z = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  double psi = 0.5;

  for (std::size_t v = 0; v < 9; ++v) {
    auto probs = mdgm.ZFullConditional(z, v, psi);
    EXPECT_EQ(probs.size(), 3);
    double sum = 0.0;
    for (double p : probs) sum += p;
    EXPECT_NEAR(sum, 1.0, 1e-10) << "vertex=" << v;
    for (double p : probs) {
      EXPECT_GE(p, 0.0) << "vertex=" << v;
    }
  }
}

TEST(MixtureDirectedGraphicalModels, FourColorLogLikelihoodFinite) {
  auto nug = MakeGrid3x3();
  MixtureDirectedGraphicalModels mdgm(nug, DagType::kSpanningTree, 4);

  std::vector<int> z = {0, 1, 2, 3, 0, 1, 2, 3, 0};
  double ll = mdgm.LogLikelihood(z, 0.5);
  EXPECT_TRUE(std::isfinite(ll));
  EXPECT_LT(ll, 0.0);
}

TEST(EmissionBernoulli, ThreeColorUpdateMaintainsOrdering) {
  // Build observations: 4 vertices, 3 colors
  Observations y;
  // vertex 0: {1, 0, 1}, vertex 1: {0, 0}, vertex 2: {1, 1, 1}, vertex 3: {0, 1}
  y.data = {1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0};
  y.ptr = {0, 3, 5, 8, 10};

  std::vector<int> z = {0, 0, 2, 1};
  std::vector<double> theta = {0.2, 0.5, 0.8};
  std::vector<double> prior = {1.0, 1.0};
  RNG rng(42);

  for (int iter = 0; iter < 100; ++iter) {
    theta = UpdateEmissionParams(y, z, theta, prior, 3, FamilyType::kBernoulli, rng);
    EXPECT_EQ(theta.size(), 3);
    // Verify ordering constraint
    EXPECT_LT(theta[0], theta[1]) << "iter=" << iter;
    EXPECT_LT(theta[1], theta[2]) << "iter=" << iter;
    // All in [0, 1]
    for (double e : theta) {
      EXPECT_GE(e, 0.0) << "iter=" << iter;
      EXPECT_LE(e, 1.0) << "iter=" << iter;
    }
  }
}

TEST(EmissionBernoulli, FourColorUpdateMaintainsOrdering) {
  Observations y;
  y.data = {1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0};
  y.ptr = {0, 3, 5, 8, 10, 12};

  std::vector<int> z = {0, 1, 3, 2, 1};
  std::vector<double> theta = {0.1, 0.3, 0.6, 0.9};
  std::vector<double> prior = {1.0, 1.0};
  RNG rng(123);

  for (int iter = 0; iter < 100; ++iter) {
    theta = UpdateEmissionParams(y, z, theta, prior, 4, FamilyType::kBernoulli, rng);
    EXPECT_EQ(theta.size(), 4);
    for (std::size_t k = 0; k + 1 < theta.size(); ++k) {
      EXPECT_LT(theta[k], theta[k + 1]) << "iter=" << iter << " k=" << k;
    }
  }
}

TEST(EmissionBernoulli, ThreeColorLikelihoodCorrect) {
  // Single vertex with observations {1, 0}
  std::vector<double> obs = {1.0, 0.0};
  std::vector<double> theta = {0.2, 0.5, 0.8};

  auto lik = EmissionLikelihood(obs, theta, 3, FamilyType::kBernoulli);
  EXPECT_EQ(lik.size(), 3);
  // p(y={1,0} | k=0) = 0.2 * 0.8 = 0.16
  EXPECT_NEAR(lik[0], 0.2 * 0.8, 1e-10);
  // p(y={1,0} | k=1) = 0.5 * 0.5 = 0.25
  EXPECT_NEAR(lik[1], 0.5 * 0.5, 1e-10);
  // p(y={1,0} | k=2) = 0.8 * 0.2 = 0.16
  EXPECT_NEAR(lik[2], 0.8 * 0.2, 1e-10);
}

TEST(HierarchicalModel, ThreeColorZFullConditionalWithEmission) {
  // Chain graph: 0-1-2
  std::vector<std::size_t> row = {0, 1, 1, 2};
  std::vector<std::size_t> col = {1, 0, 2, 1};
  NaturalUndirectedGraph nug(GraphCOO(3, row, col));

  auto spatial = std::make_unique<MixtureDirectedGraphicalModels>(
      nug, DagType::kSpanningTree, 3);
  Model model(std::move(spatial), FamilyType::kBernoulli);

  Observations y;
  y.data = {1.0, 1.0, 0.0, 0.0, 1.0, 0.0};
  y.ptr = {0, 2, 4, 6};

  std::vector<int> z = {0, 1, 2};
  std::vector<double> theta = {0.2, 0.5, 0.8};

  auto probs = model.ZFullConditional(z, 1, 1.0, y, theta);
  EXPECT_EQ(probs.size(), 3);
  double sum = 0.0;
  for (double p : probs) sum += p;
  EXPECT_NEAR(sum, 1.0, 1e-10);
  for (double p : probs) {
    EXPECT_GE(p, 0.0);
  }
}

// --- Gaussian emission tests ---

TEST(EmissionGaussian, LikelihoodCorrect) {
  // Single vertex with observation y=5
  std::vector<double> obs = {5.0};
  // theta: [mu_0, mu_1, sigma2_0, sigma2_1]
  std::vector<double> theta = {3.0, 7.0, 1.0, 4.0};

  auto lik = EmissionLikelihood(obs, theta, 2, FamilyType::kGaussian);
  EXPECT_EQ(lik.size(), 2);
  // N(5 | 3, sigma2=1) = exp(-2) / sqrt(2*pi)
  double expected_0 = std::exp(-0.5 * 4.0) / std::sqrt(2.0 * M_PI);
  // N(5 | 7, sigma2=4) = exp(-0.5 * 1) / (2 * sqrt(2*pi))
  double expected_1 = std::exp(-0.5 * 1.0) / (2.0 * std::sqrt(2.0 * M_PI));
  EXPECT_NEAR(lik[0], expected_0, 1e-10);
  EXPECT_NEAR(lik[1], expected_1, 1e-10);
}

TEST(EmissionGaussian, UpdateMaintainsMuOrdering) {
  Observations y;
  y.data = {1.0, 2.0, 3.0, 8.0, 9.0, 10.0};
  y.ptr = {0, 3, 6};

  std::vector<int> z = {0, 1};
  // theta: [mu_0, mu_1, sigma2_0, sigma2_1]
  std::vector<double> theta = {2.0, 9.0, 1.0, 1.0};
  std::vector<double> prior = {0.0, 10000.0, 2.0, 1.0};  // mu_0, sigma2_0, alpha_0, beta_0
  RNG rng(42);

  for (int iter = 0; iter < 50; ++iter) {
    theta = UpdateEmissionParams(y, z, theta, prior, 2, FamilyType::kGaussian, rng);
    EXPECT_EQ(theta.size(), 4);
    // mu ordering: mu[0] < mu[1]
    EXPECT_LT(theta[0], theta[1]) << "iter=" << iter;
    // sigma2 must be positive
    EXPECT_GT(theta[2], 0.0) << "iter=" << iter;
    EXPECT_GT(theta[3], 0.0) << "iter=" << iter;
  }
}

// --- Poisson emission tests ---

TEST(EmissionPoisson, LikelihoodCorrect) {
  std::vector<double> obs = {3.0};
  std::vector<double> theta = {1.0, 5.0};

  auto lik = EmissionLikelihood(obs, theta, 2, FamilyType::kPoisson);
  EXPECT_EQ(lik.size(), 2);
  // Pois(3 | 1) = e^{-1} * 1^3 / 3! = e^{-1}/6
  double expected_0 = std::exp(-1.0) / 6.0;
  // Pois(3 | 5) = e^{-5} * 5^3 / 3! = 125 * e^{-5} / 6
  double expected_1 = 125.0 * std::exp(-5.0) / 6.0;
  EXPECT_NEAR(lik[0], expected_0, 1e-10);
  EXPECT_NEAR(lik[1], expected_1, 1e-10);
}

TEST(EmissionPoisson, UpdateMaintainsOrdering) {
  Observations y;
  y.data = {0.0, 1.0, 0.0, 5.0, 6.0, 4.0};
  y.ptr = {0, 3, 6};

  std::vector<int> z = {0, 1};
  std::vector<double> theta = {1.0, 5.0};
  std::vector<double> prior = {1.0, 0.1};  // alpha_0, beta_0
  RNG rng(42);

  for (int iter = 0; iter < 50; ++iter) {
    theta = UpdateEmissionParams(y, z, theta, prior, 2, FamilyType::kPoisson, rng);
    EXPECT_EQ(theta.size(), 2);
    EXPECT_LT(theta[0], theta[1]) << "iter=" << iter;
    EXPECT_GT(theta[0], 0.0) << "iter=" << iter;
    EXPECT_GT(theta[1], 0.0) << "iter=" << iter;
  }
}

TEST(EmissionPoisson, ThreeColorUpdateMaintainsOrdering) {
  Observations y;
  y.data = {0.0, 1.0, 3.0, 4.0, 8.0, 9.0};
  y.ptr = {0, 2, 4, 6};

  std::vector<int> z = {0, 1, 2};
  std::vector<double> theta = {0.5, 3.5, 8.5};
  std::vector<double> prior = {1.0, 0.1};
  RNG rng(42);

  for (int iter = 0; iter < 50; ++iter) {
    theta = UpdateEmissionParams(y, z, theta, prior, 3, FamilyType::kPoisson, rng);
    EXPECT_EQ(theta.size(), 3);
    for (std::size_t k = 0; k + 1 < theta.size(); ++k) {
      EXPECT_LT(theta[k], theta[k + 1]) << "iter=" << iter << " k=" << k;
    }
  }
}

}  // namespace mdgm
