#include <cmath>
#include <cstddef>
#include <mdgm/markov_random_field.hpp>
#include <mdgm/mcmc.hpp>
#include <mdgm/natural_undirected_graph.hpp>
#include <mdgm/rng.hpp>
#include <span>
#include <vector>

namespace mdgm {

MarkovRandomField::MarkovRandomField(NaturalUndirectedGraph nug,
                                     PsiMethod method,
                                     std::size_t n_colors,
                                     std::size_t n_aux_sweeps)
    : nug_(std::move(nug)),
      method_(method),
      ncolors_(n_colors),
      n_aux_sweeps_(n_aux_sweeps) {}

std::vector<double> MarkovRandomField::ZFullConditional(
    std::span<const int> z, std::size_t vertex, double psi) const {
  std::vector<double> probs(ncolors_, 0.0);
  auto nbrs = nug_.neighbors(vertex);

  for (std::size_t k = 0; k < ncolors_; ++k) {
    double count = 0.0;
    for (std::size_t nb : nbrs) {
      count += static_cast<double>(static_cast<int>(k) == z[nb]);
    }
    probs[k] = std::exp(psi * count);
  }

  // Normalize
  double sum = 0.0;
  for (double p : probs) sum += p;
  if (sum > 0.0) {
    for (double& p : probs) p /= sum;
  }

  return probs;
}

double MarkovRandomField::LogLikelihood(
    std::span<const int> z, double psi) const {
  // Pseudo-log-likelihood: sum_i log p(z_i | z_{-i}, psi)
  double pll = 0.0;
  const std::size_t n = nug_.nvertices();

  for (std::size_t i = 0; i < n; ++i) {
    auto nbrs = nug_.neighbors(i);

    // Count same-color neighbors for vertex i's actual color
    double n_same = 0.0;
    for (std::size_t nb : nbrs) {
      n_same += static_cast<double>(z[i] == z[nb]);
    }

    // Log-sum-exp for normalizing constant
    double max_log = -1e300;
    std::vector<double> log_terms(ncolors_);
    for (std::size_t k = 0; k < ncolors_; ++k) {
      double count = 0.0;
      for (std::size_t nb : nbrs) {
        count += static_cast<double>(static_cast<int>(k) == z[nb]);
      }
      log_terms[k] = psi * count;
      if (log_terms[k] > max_log) max_log = log_terms[k];
    }
    double lse = 0.0;
    for (std::size_t k = 0; k < ncolors_; ++k) {
      lse += std::exp(log_terms[k] - max_log);
    }

    pll += psi * n_same - max_log - std::log(lse);
  }

  return pll;
}

void MarkovRandomField::UpdateGraph(
    std::span<const int> /*z*/, double /*psi*/, RNG& /*rng*/) {
  // No-op: MRF has a fixed graph structure
}

void MarkovRandomField::StoreSample(
    std::vector<std::size_t>& dag_data,
    std::size_t iteration, std::size_t n) const {
  // Fill with SIZE_MAX (becomes NA in R) — no DAG for MRF
  for (std::size_t i = 0; i < n; ++i) {
    dag_data[i + iteration * n] = SIZE_MAX;
  }
}

double MarkovRandomField::UpdatePsi(
    std::span<const int> z, double psi,
    double psi_tune, std::size_t& accepted, RNG& rng) {
  if (method_ == PsiMethod::kPseudoLikelihood) {
    // Use default MH with pseudo-log-likelihood
    return SpatialRandomField::UpdatePsi(z, psi, psi_tune, accepted, rng);
  }

  // Exchange algorithm
  double proposal = psi + rng.normal(0.0, psi_tune);
  if (proposal <= 0.0) return psi;

  // Sample auxiliary field from the MRF at the proposed psi
  std::vector<int> z_aux = (ncolors_ == 2) ? CftpSample(proposal, rng)
                                           : GibbsSample(proposal, rng);

  // Sufficient statistics
  double S_z = SufficientStatistic(z);
  double S_aux = SufficientStatistic(z_aux);

  // Exchange ratio (partition functions cancel)
  double log_ratio = S_z * proposal + S_aux * psi
                   - S_z * psi - S_aux * proposal
                   + LogPriorPsi(proposal) - LogPriorPsi(psi);

  if (std::log(rng.uniform()) < log_ratio) {
    ++accepted;
    return proposal;
  }
  return psi;
}

double MarkovRandomField::SufficientStatistic(
    std::span<const int> z) const {
  // Count same-color neighbor pairs (each undirected edge counted once)
  double count = 0.0;
  const std::size_t n = nug_.nvertices();
  for (std::size_t i = 0; i < n; ++i) {
    auto nbrs = nug_.neighbors(i);
    for (std::size_t nb : nbrs) {
      if (nb > i) {  // count each edge once
        count += static_cast<double>(z[i] == z[nb]);
      }
    }
  }
  return count;
}

std::vector<int> MarkovRandomField::GibbsSample(
    double psi, RNG& rng) const {
  const std::size_t n = nug_.nvertices();

  // Initialize with random coloring
  std::vector<int> z(n);
  for (std::size_t i = 0; i < n; ++i) {
    z[i] = static_cast<int>(rng.uniform<std::size_t>(0, ncolors_ - 1));
  }

  // Gibbs sweeps
  for (std::size_t sweep = 0; sweep < n_aux_sweeps_; ++sweep) {
    auto perm = rng.permutation(n);
    for (std::size_t idx = 0; idx < n; ++idx) {
      std::size_t i = perm[idx];
      auto probs = ZFullConditional(z, i, psi);
      z[i] = static_cast<int>(rng.discrete(probs));
    }
  }

  return z;
}

std::vector<int> MarkovRandomField::Sample(
    double psi, RNG& rng) const {
  return (ncolors_ == 2) ? CftpSample(psi, rng) : GibbsSample(psi, rng);
}

void MarkovRandomField::CftpSweep(
    std::vector<int>& config, double psi,
    const std::vector<double>& uniforms,
    std::size_t col_offset) const {
  const std::size_t n = nug_.nvertices();
  for (std::size_t i = 0; i < n; ++i) {
    auto nbrs = nug_.neighbors(i);
    double n0 = 0.0;
    double n1 = 0.0;
    for (std::size_t nb : nbrs) {
      if (config[nb] == 0) n0 += 1.0;
      else n1 += 1.0;
    }
    double w0 = std::exp(psi * n0);
    double w1 = std::exp(psi * n1);
    double threshold = w0 / (w0 + w1);
    config[i] = (uniforms[col_offset * n + i] < threshold) ? 0 : 1;
  }
}

std::vector<int> MarkovRandomField::CftpSample(
    double psi, RNG& rng) const {
  const std::size_t n = nug_.nvertices();
  std::size_t block_size = n_aux_sweeps_;
  constexpr std::size_t kMaxDoublings = 30;

  // Shared uniforms: flat vector of size n * total_sweeps
  // uniforms[sweep * n + site]
  std::vector<double> uniforms(n * block_size);
  for (auto& u : uniforms) u = rng.uniform();

  for (std::size_t doubling = 0; doubling <= kMaxDoublings; ++doubling) {
    std::size_t total_sweeps = uniforms.size() / n;

    // Initialize extremal states
    std::vector<int> lo(n, 0);
    std::vector<int> hi(n, 1);

    // Replay all sweeps with shared uniforms
    for (std::size_t s = 0; s < total_sweeps; ++s) {
      CftpSweep(lo, psi, uniforms, s);
      CftpSweep(hi, psi, uniforms, s);
    }

    // Check coalescence
    if (lo == hi) return lo;

    // Double by prepending new uniforms
    std::vector<double> new_uniforms(n * total_sweeps);
    for (auto& u : new_uniforms) u = rng.uniform();
    new_uniforms.insert(new_uniforms.end(), uniforms.begin(), uniforms.end());
    uniforms = std::move(new_uniforms);
  }

  // Safety fallback: return last lower chain (should not happen in practice)
  std::size_t total_sweeps = uniforms.size() / n;
  std::vector<int> lo(n, 0);
  for (std::size_t s = 0; s < total_sweeps; ++s) {
    CftpSweep(lo, psi, uniforms, s);
  }
  return lo;
}

std::size_t MarkovRandomField::nvertices() const {
  return nug_.nvertices();
}

std::size_t MarkovRandomField::ncolors() const {
  return ncolors_;
}

}  // namespace mdgm
