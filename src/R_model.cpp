#include <cpp11.hpp>
#include <cstddef>
#include <mdgm/directed_acyclic_graph.hpp>
#include <mdgm/markov_random_field.hpp>
#include <mdgm/mcmc.hpp>
#include <mdgm/mixture_directed_graphical_models.hpp>
#include <mdgm/model.hpp>
#include <mdgm/natural_undirected_graph.hpp>
#include <mdgm/observations.hpp>
#include <mdgm/rng.hpp>
#include <mdgm/spatial_random_field.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

inline mdgm::DagType ParseDagType(const std::string& type) {
  if (type == "spanning_tree") return mdgm::DagType::kSpanningTree;
  if (type == "acyclic_orientation") return mdgm::DagType::kAcyclicOrientation;
  if (type == "rooted") return mdgm::DagType::kRooted;
  throw std::invalid_argument("Unknown DAG type: " + type);
}

inline mdgm::PsiMethod ParsePsiMethod(const std::string& method) {
  if (method == "exchange") return mdgm::PsiMethod::kExchange;
  if (method == "pseudo_likelihood") return mdgm::PsiMethod::kPseudoLikelihood;
  throw std::invalid_argument("Unknown MRF method: " + method);
}

inline mdgm::FamilyType ParseFamilyType(const std::string& type) {
  if (type == "bernoulli") return mdgm::FamilyType::kBernoulli;
  if (type == "gaussian") return mdgm::FamilyType::kGaussian;
  if (type == "poisson") return mdgm::FamilyType::kPoisson;
  if (type == "categorical") return mdgm::FamilyType::kCategorical;
  throw std::invalid_argument("Unknown emission family: " + type);
}

// Build Observations from R flat vectors (data + ptr).
// ptr is 0-indexed offsets of length nvertices+1; data is flat numeric vector.
inline mdgm::Observations BuildObservations(
    const cpp11::doubles& obs_data, const cpp11::integers& obs_ptr) {
  mdgm::Observations y;
  y.data.assign(obs_data.begin(), obs_data.end());
  y.ptr.resize(static_cast<std::size_t>(obs_ptr.size()));
  for (std::size_t i = 0; i < y.ptr.size(); ++i) {
    y.ptr[i] = static_cast<std::size_t>(obs_ptr[static_cast<R_xlen_t>(i)]);
  }
  return y;
}

}  // namespace

// --- Model creation ---

[[cpp11::register]]
cpp11::external_pointer<mdgm::Model> model_create_standalone_cpp(
    cpp11::external_pointer<mdgm::NaturalUndirectedGraph> nug,
    cpp11::strings dag_type, int n_colors) {
  auto dt = ParseDagType(static_cast<std::string>(dag_type[0]));
  auto spatial = std::make_unique<mdgm::MixtureDirectedGraphicalModels>(
      *nug, dt, static_cast<std::size_t>(n_colors));
  auto model = std::make_unique<mdgm::Model>(std::move(spatial));
  return cpp11::external_pointer<mdgm::Model>(model.release());
}

[[cpp11::register]]
cpp11::external_pointer<mdgm::Model> model_create_hierarchical_cpp(
    cpp11::external_pointer<mdgm::NaturalUndirectedGraph> nug,
    cpp11::strings dag_type, int n_colors, cpp11::strings emission) {
  auto dt = ParseDagType(static_cast<std::string>(dag_type[0]));
  auto ft = ParseFamilyType(static_cast<std::string>(emission[0]));
  auto spatial = std::make_unique<mdgm::MixtureDirectedGraphicalModels>(
      *nug, dt, static_cast<std::size_t>(n_colors));
  auto model = std::make_unique<mdgm::Model>(std::move(spatial), ft);
  return cpp11::external_pointer<mdgm::Model>(model.release());
}

// --- MRF model creation ---

[[cpp11::register]]
cpp11::external_pointer<mdgm::Model> model_create_mrf_standalone_cpp(
    cpp11::external_pointer<mdgm::NaturalUndirectedGraph> nug,
    cpp11::strings method, int n_colors, int n_aux_sweeps) {
  auto pm = ParsePsiMethod(static_cast<std::string>(method[0]));
  auto spatial = std::make_unique<mdgm::MarkovRandomField>(
      *nug, pm, static_cast<std::size_t>(n_colors),
      static_cast<std::size_t>(n_aux_sweeps));
  auto model = std::make_unique<mdgm::Model>(std::move(spatial));
  return cpp11::external_pointer<mdgm::Model>(model.release());
}

[[cpp11::register]]
cpp11::external_pointer<mdgm::Model> model_create_mrf_hierarchical_cpp(
    cpp11::external_pointer<mdgm::NaturalUndirectedGraph> nug,
    cpp11::strings method, int n_colors, cpp11::strings emission,
    int n_aux_sweeps) {
  auto pm = ParsePsiMethod(static_cast<std::string>(method[0]));
  auto ft = ParseFamilyType(static_cast<std::string>(emission[0]));
  auto spatial = std::make_unique<mdgm::MarkovRandomField>(
      *nug, pm, static_cast<std::size_t>(n_colors),
      static_cast<std::size_t>(n_aux_sweeps));
  auto model = std::make_unique<mdgm::Model>(std::move(spatial), ft);
  return cpp11::external_pointer<mdgm::Model>(model.release());
}

[[cpp11::register]]
bool model_has_emission_cpp(cpp11::external_pointer<mdgm::Model> model) {
  return model->has_emission();
}

[[cpp11::register]]
int model_nvertices_cpp(cpp11::external_pointer<mdgm::Model> model) {
  return static_cast<int>(model->nvertices());
}

[[cpp11::register]]
int model_ncolors_cpp(cpp11::external_pointer<mdgm::Model> model) {
  return static_cast<int>(model->ncolors());
}

[[cpp11::register]]
cpp11::sexp model_emission_type_cpp(cpp11::external_pointer<mdgm::Model> model) {
  if (!model->has_emission()) {
    return R_NilValue;
  }
  // Use the FamilyType enum to return a string
  auto ft = model->emission_type();
  switch (ft) {
    case mdgm::FamilyType::kBernoulli: return cpp11::as_sexp("bernoulli");
    case mdgm::FamilyType::kGaussian: return cpp11::as_sexp("gaussian");
    case mdgm::FamilyType::kPoisson: return cpp11::as_sexp("poisson");
    case mdgm::FamilyType::kCategorical: return cpp11::as_sexp("categorical");
    default: return R_NilValue;
  }
}

// --- MCMC ---

[[cpp11::register]]
cpp11::writable::list run_mcmc_cpp(
    cpp11::external_pointer<mdgm::Model> model,
    const cpp11::doubles& obs_data, const cpp11::integers& obs_ptr,
    const cpp11::integers& z_init, double psi_init,
    const cpp11::doubles& theta_init,
    int n_iterations, double psi_tune,
    const cpp11::doubles& emission_prior_params,
    cpp11::external_pointer<mdgm::RNG> rng) {
  using namespace cpp11::literals;

  // Build Observations
  mdgm::Observations y = BuildObservations(obs_data, obs_ptr);

  // Build z_init (convert from 0-indexed R integers)
  std::vector<int> z0(z_init.begin(), z_init.end());

  // Build theta_init
  std::vector<double> theta0(theta_init.begin(), theta_init.end());

  // Build config
  mdgm::McmcConfig config;
  config.n_iterations = static_cast<std::size_t>(n_iterations);
  config.psi_tune = psi_tune;
  config.emission_prior_params.assign(
      emission_prior_params.begin(), emission_prior_params.end());

  // Run
  mdgm::McmcSamples samples =
      mdgm::RunMcmc(*model, y, z0, psi_init, theta0, config, *rng);

  // Convert to R objects
  const std::size_t n = samples.n_vertices;
  const std::size_t J = samples.n_iterations;
  const std::size_t nc = samples.n_colors;
  const std::size_t n_theta = samples.n_theta;

  // z: integer matrix n x J (column-major, already stored that way)
  cpp11::writable::integers z_r(static_cast<R_xlen_t>(n * J));
  for (std::size_t i = 0; i < n * J; ++i) {
    z_r[static_cast<R_xlen_t>(i)] = samples.z[i];
  }
  z_r.attr("dim") = cpp11::writable::integers({
      static_cast<int>(n), static_cast<int>(J)});

  // psi: numeric vector length J
  cpp11::writable::doubles psi_r(static_cast<R_xlen_t>(J));
  for (std::size_t j = 0; j < J; ++j) {
    psi_r[static_cast<R_xlen_t>(j)] = samples.psi[j];
  }

  // theta: numeric matrix n_theta x J (only if emission model present)
  cpp11::sexp theta_r;
  if (samples.theta.empty()) {
    theta_r = R_NilValue;
  } else {
    cpp11::writable::doubles theta_tmp(static_cast<R_xlen_t>(n_theta * J));
    for (std::size_t i = 0; i < n_theta * J; ++i) {
      theta_tmp[static_cast<R_xlen_t>(i)] = samples.theta[i];
    }
    theta_tmp.attr("dim") = cpp11::writable::integers({
        static_cast<int>(n_theta), static_cast<int>(J)});
    theta_r = theta_tmp;
  }

  // dag_data: integer matrix n x J
  cpp11::writable::integers dag_r(static_cast<R_xlen_t>(n * J));
  for (std::size_t i = 0; i < n * J; ++i) {
    // Convert SIZE_MAX (root marker) to NA, otherwise 1-indexed
    if (samples.dag_data[i] == std::numeric_limits<std::size_t>::max()) {
      dag_r[static_cast<R_xlen_t>(i)] = NA_INTEGER;
    } else {
      dag_r[static_cast<R_xlen_t>(i)] =
          static_cast<int>(samples.dag_data[i] + 1);
    }
  }
  dag_r.attr("dim") = cpp11::writable::integers({
      static_cast<int>(n), static_cast<int>(J)});

  return cpp11::writable::list({
      "z"_nm = z_r,
      "psi"_nm = psi_r,
      "theta"_nm = theta_r,
      "dag"_nm = dag_r,
      "psi_accepted"_nm = static_cast<int>(samples.psi_accepted),
      "graph_accepted"_nm = static_cast<int>(samples.graph_accepted),
      "n_iterations"_nm = static_cast<int>(J),
      "n_vertices"_nm = static_cast<int>(n),
      "n_colors"_nm = static_cast<int>(nc)
  });
}
