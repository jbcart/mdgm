#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <mdgm/graph_storage.hpp>
#include <mdgm/markov_random_field.hpp>
#include <mdgm/model.hpp>
#include <mdgm/natural_undirected_graph.hpp>
#include <mdgm/rng.hpp>
#include <memory>
#include <numeric>
#include <vector>

namespace mdgm {

// Helper: chain graph 0-1-2
static NaturalUndirectedGraph MakeChain3() {
  std::vector<std::size_t> row = {0, 1, 1, 2};
  std::vector<std::size_t> col = {1, 0, 2, 1};
  return NaturalUndirectedGraph(GraphCOO(3, row, col));
}

// Helper: triangle graph 0-1-2-0
static NaturalUndirectedGraph MakeTriangleMRF() {
  std::vector<std::size_t> row = {0, 0, 1, 1, 2, 2};
  std::vector<std::size_t> col = {1, 2, 0, 2, 0, 1};
  return NaturalUndirectedGraph(GraphCOO(3, row, col));
}

// Helper: 3x3 grid
static NaturalUndirectedGraph MakeGrid3x3MRF() {
  return GenerateRegularGraph({3, 3}, 1);
}

TEST(MarkovRandomField, Construction) {
  auto nug = MakeTriangleMRF();
  EXPECT_NO_THROW({
    MarkovRandomField mrf(nug, PsiMethod::kExchange, 2, 100);
  });
  MarkovRandomField mrf(nug, PsiMethod::kPseudoLikelihood, 3, 50);
  EXPECT_EQ(mrf.nvertices(), 3);
  EXPECT_EQ(mrf.ncolors(), 3);
}

TEST(MarkovRandomField, ZFullConditionalSumsToOne) {
  auto nug = MakeGrid3x3MRF();
  MarkovRandomField mrf(nug, PsiMethod::kExchange, 2);

  std::vector<int> z = {0, 1, 0, 1, 0, 1, 0, 1, 0};
  double psi = 0.5;

  for (std::size_t v = 0; v < 9; ++v) {
    auto probs = mrf.ZFullConditional(z, v, psi);
    EXPECT_EQ(probs.size(), 2);
    double sum = 0.0;
    for (double p : probs) sum += p;
    EXPECT_NEAR(sum, 1.0, 1e-10) << "vertex=" << v;
    for (double p : probs) {
      EXPECT_GE(p, 0.0) << "vertex=" << v;
    }
  }
}

TEST(MarkovRandomField, ZFullConditionalChain) {
  auto nug = MakeChain3();
  MarkovRandomField mrf(nug, PsiMethod::kExchange, 2);

  // z = {0, ?, 1}, vertex 1 has neighbors 0 and 2
  // p(z_1=0) prop to exp(psi * 1) = exp(1) [matches z[0]=0]
  // p(z_1=1) prop to exp(psi * 1) = exp(1) [matches z[2]=1]
  std::vector<int> z = {0, 0, 1};
  double psi = 1.0;
  auto probs = mrf.ZFullConditional(z, 1, psi);
  EXPECT_NEAR(probs[0], 0.5, 1e-10);
  EXPECT_NEAR(probs[1], 0.5, 1e-10);

  // z = {0, ?, 0}: both neighbors are 0
  // p(z_1=0) prop to exp(2*psi) = exp(2)
  // p(z_1=1) prop to exp(0) = 1
  z = {0, 0, 0};
  probs = mrf.ZFullConditional(z, 1, psi);
  double e2 = std::exp(2.0);
  EXPECT_NEAR(probs[0], e2 / (e2 + 1.0), 1e-10);
  EXPECT_NEAR(probs[1], 1.0 / (e2 + 1.0), 1e-10);
}

TEST(MarkovRandomField, ZFullConditionalThreeColors) {
  auto nug = MakeGrid3x3MRF();
  MarkovRandomField mrf(nug, PsiMethod::kExchange, 3);

  std::vector<int> z = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  for (std::size_t v = 0; v < 9; ++v) {
    auto probs = mrf.ZFullConditional(z, v, 0.5);
    EXPECT_EQ(probs.size(), 3);
    double sum = 0.0;
    for (double p : probs) sum += p;
    EXPECT_NEAR(sum, 1.0, 1e-10);
  }
}

TEST(MarkovRandomField, LogLikelihoodFinite) {
  auto nug = MakeGrid3x3MRF();
  MarkovRandomField mrf(nug, PsiMethod::kPseudoLikelihood, 2);

  std::vector<int> z = {0, 1, 0, 1, 0, 1, 0, 1, 0};
  double pll = mrf.LogLikelihood(z, 0.5);
  EXPECT_TRUE(std::isfinite(pll));
  EXPECT_LT(pll, 0.0);
}

TEST(MarkovRandomField, LogLikelihoodZeroPsi) {
  // psi=0: all full conditionals uniform, pseudo-ll = n * log(1/K)
  auto nug = MakeGrid3x3MRF();
  MarkovRandomField mrf(nug, PsiMethod::kPseudoLikelihood, 2);

  std::vector<int> z = {0, 1, 0, 1, 0, 1, 0, 1, 0};
  double pll = mrf.LogLikelihood(z, 0.0);
  EXPECT_NEAR(pll, -9.0 * std::log(2.0), 1e-10);
}

TEST(MarkovRandomField, UpdateGraphIsNoOp) {
  auto nug = MakeTriangleMRF();
  MarkovRandomField mrf(nug, PsiMethod::kExchange, 2);
  RNG rng(42);

  std::vector<int> z = {0, 1, 0};
  // Should not throw or change anything
  EXPECT_NO_THROW(mrf.UpdateGraph(z, 1.0, rng));
}

TEST(MarkovRandomField, StoreSampleAllNA) {
  auto nug = MakeTriangleMRF();
  MarkovRandomField mrf(nug, PsiMethod::kExchange, 2);

  std::size_t n = 3;
  std::vector<std::size_t> dag_data(n * 2, 0);
  mrf.StoreSample(dag_data, 0, n);

  for (std::size_t i = 0; i < n; ++i) {
    EXPECT_EQ(dag_data[i], SIZE_MAX);
  }
}

TEST(MarkovRandomField, UpdatePsiExchangeRuns) {
  auto nug = MakeGrid3x3MRF();
  // Use small n_aux_sweeps for speed
  MarkovRandomField mrf(nug, PsiMethod::kExchange, 2, 10);
  RNG rng(42);

  std::vector<int> z = {0, 0, 0, 1, 1, 1, 0, 0, 0};
  double psi = 1.0;
  std::size_t accepted = 0;

  // Run several updates — should not crash
  for (int i = 0; i < 20; ++i) {
    psi = mrf.UpdatePsi(z, psi, 0.3, accepted, rng);
    EXPECT_GT(psi, 0.0);
    EXPECT_TRUE(std::isfinite(psi));
  }
}

TEST(MarkovRandomField, UpdatePsiPseudoLikelihoodRuns) {
  auto nug = MakeGrid3x3MRF();
  MarkovRandomField mrf(nug, PsiMethod::kPseudoLikelihood, 2);
  RNG rng(42);

  std::vector<int> z = {0, 0, 0, 1, 1, 1, 0, 0, 0};
  double psi = 1.0;
  std::size_t accepted = 0;

  for (int i = 0; i < 50; ++i) {
    psi = mrf.UpdatePsi(z, psi, 0.3, accepted, rng);
    EXPECT_TRUE(std::isfinite(psi));
  }
  // Should have some acceptance
  EXPECT_GT(accepted, static_cast<std::size_t>(0));
}

TEST(MarkovRandomField, GibbsSampleReturnsValidColoring) {
  auto nug = MakeGrid3x3MRF();
  // Access GibbsSample through UpdatePsi exchange method indirectly
  // Test that exchange algorithm produces valid psi updates
  MarkovRandomField mrf(nug, PsiMethod::kExchange, 3, 20);
  RNG rng(42);

  std::vector<int> z = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  double psi = 0.5;
  std::size_t accepted = 0;

  for (int i = 0; i < 10; ++i) {
    double new_psi = mrf.UpdatePsi(z, psi, 0.5, accepted, rng);
    EXPECT_TRUE(std::isfinite(new_psi));
    psi = new_psi;
  }
}

TEST(MarkovRandomField, ModelIntegrationStandalone) {
  auto nug = MakeGrid3x3MRF();
  auto spatial = std::make_unique<MarkovRandomField>(
      nug, PsiMethod::kPseudoLikelihood, 2);
  Model model(std::move(spatial));

  EXPECT_FALSE(model.has_emission());
  EXPECT_EQ(model.nvertices(), 9);
  EXPECT_EQ(model.ncolors(), 2);
}

TEST(MarkovRandomField, ModelIntegrationHierarchical) {
  auto nug = MakeGrid3x3MRF();
  auto spatial = std::make_unique<MarkovRandomField>(
      nug, PsiMethod::kExchange, 2, 10);
  Model model(std::move(spatial), FamilyType::kBernoulli);

  EXPECT_TRUE(model.has_emission());
  EXPECT_EQ(model.emission_type(), FamilyType::kBernoulli);
}

TEST(MarkovRandomField, CftpReturnsValidBinaryColoring) {
  auto nug = MakeGrid3x3MRF();
  MarkovRandomField mrf(nug, PsiMethod::kExchange, 2, 50);
  RNG rng(123);

  auto z = mrf.Sample(0.5, rng);
  EXPECT_EQ(z.size(), 9);
  for (int val : z) {
    EXPECT_GE(val, 0);
    EXPECT_LE(val, 1);
  }
}

TEST(MarkovRandomField, CftpCoalescesSmallGrid) {
  // Small grid with moderate psi should coalesce quickly
  auto nug = GenerateRegularGraph({4, 4}, 1);
  MarkovRandomField mrf(nug, PsiMethod::kExchange, 2, 100);
  RNG rng(42);

  // Should succeed (not hit safety limit) for multiple draws
  for (int rep = 0; rep < 5; ++rep) {
    auto z = mrf.Sample(0.3, rng);
    EXPECT_EQ(z.size(), 16);
    for (int val : z) {
      EXPECT_GE(val, 0);
      EXPECT_LE(val, 1);
    }
  }
}

TEST(MarkovRandomField, CftpSweepCouplingProperty) {
  // Two identical configs + shared uniforms must produce identical output
  auto nug = MakeGrid3x3MRF();
  MarkovRandomField mrf(nug, PsiMethod::kExchange, 2, 10);
  RNG rng(99);

  const std::size_t n = 9;
  std::vector<double> uniforms(n);
  for (auto& u : uniforms) u = rng.uniform();

  std::vector<int> a = {0, 1, 0, 1, 0, 1, 0, 1, 0};
  std::vector<int> b = a;  // identical

  mrf.CftpSweep(a, 0.5, uniforms, 0);
  mrf.CftpSweep(b, 0.5, uniforms, 0);

  EXPECT_EQ(a, b);
}

TEST(MarkovRandomField, SampleDispatchesByNcolors) {
  auto nug = MakeGrid3x3MRF();
  RNG rng(42);

  // Binary: uses CFTP
  MarkovRandomField mrf2(nug, PsiMethod::kExchange, 2, 50);
  auto z2 = mrf2.Sample(0.5, rng);
  EXPECT_EQ(z2.size(), 9);
  for (int val : z2) {
    EXPECT_GE(val, 0);
    EXPECT_LE(val, 1);
  }

  // 3-color: uses Gibbs
  MarkovRandomField mrf3(nug, PsiMethod::kExchange, 3, 50);
  auto z3 = mrf3.Sample(0.5, rng);
  EXPECT_EQ(z3.size(), 9);
  for (int val : z3) {
    EXPECT_GE(val, 0);
    EXPECT_LE(val, 2);
  }
}

}  // namespace mdgm
