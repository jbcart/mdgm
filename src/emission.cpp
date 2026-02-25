#include <cmath>
#include <cstddef>
#include <mdgm/emission.hpp>
#include <mdgm/observations.hpp>
#include <mdgm/rng.hpp>
#include <mdgm/spatial_random_field.hpp>
#include <span>
#include <stdexcept>
#include <vector>

namespace mdgm {

std::vector<double> EmissionLikelihood(
    std::span<const int> y_i, std::span<const double> eta,
    std::size_t ncolors, FamilyType type) {
  std::vector<double> lik(ncolors, 1.0);

  switch (type) {
    case FamilyType::kBernoulli:
      for (std::size_t k = 0; k < ncolors; ++k) {
        double p = eta[k];
        for (int y : y_i) {
          lik[k] *= (y == 1) ? p : (1.0 - p);
        }
      }
      break;
    case FamilyType::kGaussian:
      // eta layout: [mu_0, ..., mu_{K-1}, sigma_0, ..., sigma_{K-1}]
      for (std::size_t k = 0; k < ncolors; ++k) {
        double mu = eta[k];
        double sigma = eta[ncolors + k];
        for (int y : y_i) {
          double dy = static_cast<double>(y);
          double z = (dy - mu) / sigma;
          lik[k] *= std::exp(-0.5 * z * z) / (sigma * std::sqrt(2.0 * M_PI));
        }
      }
      break;
    case FamilyType::kPoisson:
      for (std::size_t k = 0; k < ncolors; ++k) {
        double lambda = eta[k];
        for (int y : y_i) {
          // Poisson PMF: lambda^y * exp(-lambda) / y!
          double log_pmf = static_cast<double>(y) * std::log(lambda) - lambda -
                           std::lgamma(static_cast<double>(y) + 1.0);
          lik[k] *= std::exp(log_pmf);
        }
      }
      break;
    default:
      throw std::invalid_argument("Unsupported emission family type");
  }

  return lik;
}

double EmissionLogLikelihood(
    const Observations& y, std::span<const int> z,
    std::span<const double> eta, FamilyType type) {
  double ll = 0.0;

  switch (type) {
    case FamilyType::kBernoulli:
      for (std::size_t i = 0; i < y.nvertices(); ++i) {
        if (y.empty(i)) continue;
        double p = eta[static_cast<std::size_t>(z[i])];
        for (int yij : y[i]) {
          ll += (yij == 1) ? std::log(p) : std::log(1.0 - p);
        }
      }
      break;
    case FamilyType::kGaussian: {
      std::size_t nc = 0;
      // Infer ncolors: eta has 2*ncolors elements
      // We need ncolors, but it's not passed. Infer from z max.
      // Actually, we need to find ncolors from eta.size() / 2
      nc = eta.size() / 2;
      for (std::size_t i = 0; i < y.nvertices(); ++i) {
        if (y.empty(i)) continue;
        std::size_t k = static_cast<std::size_t>(z[i]);
        double mu = eta[k];
        double sigma = eta[nc + k];
        for (int yij : y[i]) {
          double dy = static_cast<double>(yij);
          double zv = (dy - mu) / sigma;
          ll += -0.5 * zv * zv - std::log(sigma) - 0.5 * std::log(2.0 * M_PI);
        }
      }
      break;
    }
    case FamilyType::kPoisson:
      for (std::size_t i = 0; i < y.nvertices(); ++i) {
        if (y.empty(i)) continue;
        double lambda = eta[static_cast<std::size_t>(z[i])];
        for (int yij : y[i]) {
          ll += static_cast<double>(yij) * std::log(lambda) - lambda -
                std::lgamma(static_cast<double>(yij) + 1.0);
        }
      }
      break;
    default:
      throw std::invalid_argument("Unsupported emission family type");
  }

  return ll;
}

std::vector<double> UpdateEmissionParams(
    const Observations& y, std::span<const int> z,
    std::span<const double> eta, std::span<const double> prior_params,
    std::size_t ncolors, FamilyType type, RNG& rng) {
  std::vector<double> new_eta(eta.begin(), eta.end());

  switch (type) {
    case FamilyType::kBernoulli: {
      double a = prior_params[0];
      double b = prior_params[1];

      std::vector<double> successes(ncolors, 0.0);
      std::vector<double> trials(ncolors, 0.0);
      for (std::size_t i = 0; i < y.nvertices(); ++i) {
        std::size_t k = static_cast<std::size_t>(z[i]);
        for (int yij : y[i]) {
          trials[k] += 1.0;
          successes[k] += static_cast<double>(yij);
        }
      }

      // Sequential truncated Beta: eta[0] < eta[1] < ... < eta[K-1]
      double lower = 0.0;
      for (std::size_t k = 0; k < ncolors; ++k) {
        double upper = (k + 1 < ncolors) ? new_eta[k + 1] : 1.0;
        new_eta[k] = rng.truncated_beta(
            a + successes[k], b + trials[k] - successes[k],
            lower, upper);
        lower = new_eta[k];
      }
      break;
    }
    case FamilyType::kGaussian: {
      // Prior params: {mu_0, kappa_0, alpha_0, beta_0}
      // Normal-InverseGamma conjugate update with ordering on mu
      double mu_0 = prior_params[0];
      double kappa_0 = prior_params[1];
      double alpha_0 = prior_params[2];
      double beta_0 = prior_params[3];

      // eta layout: [mu_0, ..., mu_{K-1}, sigma_0, ..., sigma_{K-1}]
      // Sufficient statistics per color
      std::vector<double> sum_y(ncolors, 0.0);
      std::vector<double> sum_y2(ncolors, 0.0);
      std::vector<double> counts(ncolors, 0.0);

      for (std::size_t i = 0; i < y.nvertices(); ++i) {
        std::size_t k = static_cast<std::size_t>(z[i]);
        for (int yij : y[i]) {
          double dy = static_cast<double>(yij);
          sum_y[k] += dy;
          sum_y2[k] += dy * dy;
          counts[k] += 1.0;
        }
      }

      // Update sigma first (Inverse-Gamma), then mu (truncated Normal)
      for (std::size_t k = 0; k < ncolors; ++k) {
        double n_k = counts[k];
        double ybar_k = (n_k > 0) ? sum_y[k] / n_k : 0.0;
        double mu_k = new_eta[k];

        // Posterior for sigma^2: InverseGamma(alpha_n, beta_n)
        double alpha_n = alpha_0 + n_k / 2.0;
        double ss = sum_y2[k] - 2.0 * mu_k * sum_y[k] + n_k * mu_k * mu_k;
        double beta_n = beta_0 + 0.5 * ss +
                        kappa_0 * n_k * (ybar_k - mu_0) * (ybar_k - mu_0) /
                            (2.0 * (kappa_0 + n_k));
        double sigma2 = rng.inverse_gamma(alpha_n, beta_n);
        new_eta[ncolors + k] = std::sqrt(sigma2);

        // Posterior for mu: Normal(mu_n, sigma2 / kappa_n) with ordering
        double kappa_n = kappa_0 + n_k;
        double mu_n = (kappa_0 * mu_0 + n_k * ybar_k) / kappa_n;
        double mu_sd = std::sqrt(sigma2 / kappa_n);

        double lo = (k > 0) ? new_eta[k - 1] : -1e10;
        double hi = (k + 1 < ncolors) ? new_eta[k + 1] : 1e10;
        new_eta[k] = rng.truncated_normal(mu_n, mu_sd, lo, hi);
      }
      break;
    }
    case FamilyType::kPoisson: {
      // Prior params: {alpha_0, beta_0} for Gamma(alpha_0, beta_0) prior
      double alpha_0 = prior_params[0];
      double beta_0 = prior_params[1];

      std::vector<double> sum_y(ncolors, 0.0);
      std::vector<double> counts(ncolors, 0.0);

      for (std::size_t i = 0; i < y.nvertices(); ++i) {
        std::size_t k = static_cast<std::size_t>(z[i]);
        for (int yij : y[i]) {
          sum_y[k] += static_cast<double>(yij);
          counts[k] += 1.0;
        }
      }

      // Sequential truncated Gamma: lambda[0] < lambda[1] < ... < lambda[K-1]
      double lower = 0.0;
      for (std::size_t k = 0; k < ncolors; ++k) {
        double alpha_n = alpha_0 + sum_y[k];
        double beta_n = beta_0 + counts[k];
        double upper = (k + 1 < ncolors) ? new_eta[k + 1] : 1e10;
        new_eta[k] = rng.truncated_gamma(alpha_n, beta_n, lower, upper);
        lower = new_eta[k];
      }
      break;
    }
    default:
      throw std::invalid_argument("Unsupported emission family type");
  }

  return new_eta;
}

}  // namespace mdgm
