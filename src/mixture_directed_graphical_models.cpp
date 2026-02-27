#include <cmath>
#include <cstddef>
#include <mdgm/directed_acyclic_graph.hpp>
#include <mdgm/mixture_directed_graphical_models.hpp>
#include <mdgm/natural_undirected_graph.hpp>
#include <mdgm/rng.hpp>
#include <span>
#include <stdexcept>
#include <vector>

namespace mdgm {

MixtureDirectedGraphicalModels::MixtureDirectedGraphicalModels(
    NaturalUndirectedGraph nug, DagType dag_type, std::size_t n_colors)
    : nug_(std::move(nug)),
      current_dag_(GraphCOO(nug_.nvertices(), {}, {}), DagType::kSpanningTree),
      dag_type_(dag_type),
      ncolors_(n_colors),
      alpha_(n_colors, 0.0) {
  // Initialize with a proper DAG using a temporary RNG
  RNG init_rng;
  switch (dag_type_) {
    case DagType::kSpanningTree:
      current_dag_ = nug_.SampleSpanningTree(init_rng);
      break;
    case DagType::kAcyclicOrientation:
      current_dag_ = nug_.SampleAcyclicOrientation(init_rng);
      break;
    case DagType::kRooted:
      // TODO: implement rooted DAG initialization
      current_dag_ = nug_.SampleSpanningTree(init_rng);
      break;
  }
}

std::vector<double> MixtureDirectedGraphicalModels::ZFullConditional(
    std::span<const int> z, std::size_t vertex, double psi) const {
  std::vector<double> probs(ncolors_, 0.0);
  auto parents = current_dag_.parents(vertex);
  auto children = current_dag_.children(vertex);

  if (dag_type_ == DagType::kSpanningTree) {
    // Spanning trees have no v-structures, so the full conditional simplifies
    // to the conditional on neighbors (parents + children) in the tree.
    // p(z_i = k | z_{-i}) proportional to exp(psi * sum_{j in nbrs} I(k == z_j) + alpha_k)
    for (std::size_t k = 0; k < ncolors_; ++k) {
      double count = 0.0;
      for (std::size_t p : parents) {
        count += static_cast<double>(static_cast<int>(k) == z[p]);
      }
      for (std::size_t c : children) {
        count += static_cast<double>(static_cast<int>(k) == z[c]);
      }
      probs[k] = std::exp(psi * count + alpha_[k]);
    }
  } else {
    // Full conditional for AO/Rooted DAGs (Equation 5 in paper):
    // p(z_i = k | z_{-i}) proportional to
    //   exp(psi * sum_{j in parents} I(k == z_j) + alpha_k)
    //   * product over children c of:
    //     1 / sum_l exp(psi * sum_{j in parents(c)} I(l == z_j, with z_i=k) + alpha_l)

    // Numerator: contribution from parents and neighbors (parents + children)
    for (std::size_t k = 0; k < ncolors_; ++k) {
      // Count matches with parents and children (as neighbors)
      double nbr_count = 0.0;
      for (std::size_t p : parents) {
        nbr_count += static_cast<double>(static_cast<int>(k) == z[p]);
      }
      for (std::size_t c : children) {
        nbr_count += static_cast<double>(static_cast<int>(k) == z[c]);
      }
      double log_numerator = psi * nbr_count + alpha_[k];

      // Denominator: product over children of normalizing constants
      // Each child c's normalizing constant depends on z_i = k
      double log_denominator = 0.0;
      for (std::size_t c : children) {
        auto c_parents = current_dag_.parents(c);
        // Compute normalizing constant for child c when z[vertex] = k
        double max_log = -1e300;
        std::vector<double> log_terms(ncolors_);
        for (std::size_t l = 0; l < ncolors_; ++l) {
          double count = 0.0;
          for (std::size_t cp : c_parents) {
            if (cp == vertex) {
              count += static_cast<double>(static_cast<int>(l) == static_cast<int>(k));
            } else {
              count += static_cast<double>(static_cast<int>(l) == z[cp]);
            }
          }
          log_terms[l] = psi * count + alpha_[l];
          if (log_terms[l] > max_log) max_log = log_terms[l];
        }
        // Log-sum-exp for numerical stability
        double sum = 0.0;
        for (std::size_t l = 0; l < ncolors_; ++l) {
          sum += std::exp(log_terms[l] - max_log);
        }
        log_denominator += max_log + std::log(sum);
      }

      probs[k] = std::exp(log_numerator - log_denominator);
    }
  }

  // Normalize
  double sum = 0.0;
  for (double p : probs) sum += p;
  if (sum > 0.0) {
    for (double& p : probs) p /= sum;
  }

  return probs;
}

double MixtureDirectedGraphicalModels::LogLikelihood(
    std::span<const int> z, double psi) const {
  // Exact DAG log-likelihood: sum of log p(z_i | z_{parents(i)}, psi) for all i
  double ll = 0.0;
  for (std::size_t i = 0; i < nug_.nvertices(); ++i) {
    auto parents = current_dag_.parents(i);
    // Compute p(z_i | z_parents)
    double count = 0.0;
    for (std::size_t p : parents) {
      count += static_cast<double>(z[i] == z[p]);
    }
    double log_num = psi * count + alpha_[static_cast<std::size_t>(z[i])];

    // Normalizing constant
    double max_log = -1e300;
    std::vector<double> log_terms(ncolors_);
    for (std::size_t k = 0; k < ncolors_; ++k) {
      double kcount = 0.0;
      for (std::size_t p : parents) {
        kcount += static_cast<double>(static_cast<int>(k) == z[p]);
      }
      log_terms[k] = psi * kcount + alpha_[k];
      if (log_terms[k] > max_log) max_log = log_terms[k];
    }
    double sum = 0.0;
    for (std::size_t k = 0; k < ncolors_; ++k) {
      sum += std::exp(log_terms[k] - max_log);
    }
    ll += log_num - max_log - std::log(sum);
  }
  return ll;
}

void MixtureDirectedGraphicalModels::UpdateGraph(
    std::span<const int> z, double psi, RNG& rng) {
  switch (dag_type_) {
    case DagType::kSpanningTree:
      // Update weights based on z and psi, then sample spanning tree
      nug_.UpdateWeights(z, psi);
      current_dag_ = nug_.SampleSpanningTree(rng);
      break;
    case DagType::kAcyclicOrientation:
      // MH step: propose new random AO, accept/reject
      {
        DirectedAcyclicGraph proposal = nug_.SampleAcyclicOrientation(rng);
        // Temporarily swap to compute proposed likelihood
        DirectedAcyclicGraph old_dag = current_dag_;
        current_dag_ = proposal;
        double log_proposed = LogLikelihood(z, psi);
        current_dag_ = old_dag;
        double log_current = LogLikelihood(z, psi);
        double log_ratio = log_proposed - log_current;
        if (std::log(rng.uniform()) < log_ratio) {
          current_dag_ = proposal;
        }
      }
      break;
    case DagType::kRooted:
      // TODO: implement rooted DAG MH update
      break;
  }
}

void MixtureDirectedGraphicalModels::StoreSample(
    std::vector<std::size_t>& dag_data,
    std::size_t iteration, std::size_t n) const {
  // Store parent vector for compact representation
  for (std::size_t i = 0; i < n; ++i) {
    auto parents = current_dag_.parents(i);
    if (parents.empty()) {
      dag_data[i + iteration * n] = SIZE_MAX;  // root/no parent sentinel
    } else {
      dag_data[i + iteration * n] = parents[0];  // first parent (tree: exactly 1)
    }
  }
}

std::size_t MixtureDirectedGraphicalModels::nvertices() const {
  return nug_.nvertices();
}

std::size_t MixtureDirectedGraphicalModels::ncolors() const {
  return ncolors_;
}

double MixtureDirectedGraphicalModels::SufficientStatistic(
    std::span<const int> z) const {
  double count = 0.0;
  const std::size_t n = nug_.nvertices();
  for (std::size_t i = 0; i < n; ++i) {
    auto nbrs = nug_.neighbors(i);
    for (std::size_t nb : nbrs) {
      if (nb > i) {
        count += static_cast<double>(z[i] == z[nb]);
      }
    }
  }
  return count;
}

const DirectedAcyclicGraph& MixtureDirectedGraphicalModels::current_dag() const {
  return current_dag_;
}

}  // namespace mdgm
