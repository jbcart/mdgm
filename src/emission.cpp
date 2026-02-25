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
    std::span<const int> y_i, std::span<const double> theta,
    std::size_t ncolors, FamilyType type) {
  std::vector<double> lik(ncolors, 1.0);

  switch (type) {
    case FamilyType::kBernoulli:
      for (std::size_t k = 0; k < ncolors; ++k) {
        double p = theta[k];
        for (int y : y_i) {
          lik[k] *= (y == 1) ? p : (1.0 - p);
        }
      }
      break;
    case FamilyType::kGaussian:
      // theta layout: [mu_0, ..., mu_{K-1}, sigma2_0, ..., sigma2_{K-1}]
      for (std::size_t k = 0; k < ncolors; ++k) {
        double mu = theta[k];
        double sigma2 = theta[ncolors + k];
        double sigma = std::sqrt(sigma2);
        for (int y : y_i) {
          double dy = static_cast<double>(y);
          double z = (dy - mu) / sigma;
          lik[k] *= std::exp(-0.5 * z * z) / (sigma * std::sqrt(2.0 * M_PI));
        }
      }
      break;
    case FamilyType::kPoisson:
      for (std::size_t k = 0; k < ncolors; ++k) {
        double lambda = theta[k];
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
    std::span<const double> theta, FamilyType type) {
  double ll = 0.0;

  switch (type) {
    case FamilyType::kBernoulli:
      for (std::size_t i = 0; i < y.nvertices(); ++i) {
        if (y.empty(i)) continue;
        double p = theta[static_cast<std::size_t>(z[i])];
        for (int yij : y[i]) {
          ll += (yij == 1) ? std::log(p) : std::log(1.0 - p);
        }
      }
      break;
    case FamilyType::kGaussian: {
      std::size_t nc = theta.size() / 2;
      for (std::size_t i = 0; i < y.nvertices(); ++i) {
        if (y.empty(i)) continue;
        std::size_t k = static_cast<std::size_t>(z[i]);
        double mu = theta[k];
        double sigma2 = theta[nc + k];
        for (int yij : y[i]) {
          double dy = static_cast<double>(yij);
          ll += -0.5 * std::log(2.0 * M_PI * sigma2) -
                0.5 * (dy - mu) * (dy - mu) / sigma2;
        }
      }
      break;
    }
    case FamilyType::kPoisson:
      for (std::size_t i = 0; i < y.nvertices(); ++i) {
        if (y.empty(i)) continue;
        double lambda = theta[static_cast<std::size_t>(z[i])];
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
    std::span<const double> theta, std::span<const double> prior_params,
    std::size_t ncolors, FamilyType type, RNG& rng) {
  std::vector<double> new_theta(theta.begin(), theta.end());

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

      // Sequential truncated Beta: p[0] < p[1] < ... < p[K-1]
      double lower = 0.0;
      for (std::size_t k = 0; k < ncolors; ++k) {
        double upper = (k + 1 < ncolors) ? new_theta[k + 1] : 1.0;
        new_theta[k] = rng.truncated_beta(
            a + successes[k], b + trials[k] - successes[k],
            lower, upper);
        lower = new_theta[k];
      }
      break;
    }
    case FamilyType::kGaussian: {
      // Prior params: {mu_0, sigma2_0, alpha_0, beta_0}
      // Independent Normal and InverseGamma priors with ordering on mu
      // mu_k ~ N(mu_0, sigma2_0), sigma_k^2 ~ IG(alpha_0, beta_0)
      double mu_0 = prior_params[0];
      double sigma2_0 = prior_params[1];  // prior variance for mu
      double alpha_0 = prior_params[2];
      double beta_0 = prior_params[3];

      // theta layout: [mu_0, ..., mu_{K-1}, sigma2_0, ..., sigma2_{K-1}]
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
        double mu_k = new_theta[k];

        // Posterior for sigma_k^2: IG(alpha_0 + n_k/2, beta_0 + SS_k/2)
        double alpha_n = alpha_0 + n_k / 2.0;
        double ss = sum_y2[k] - 2.0 * mu_k * sum_y[k] + n_k * mu_k * mu_k;
        double beta_n = beta_0 + 0.5 * ss;
        double sigma2 = rng.inverse_gamma(alpha_n, beta_n);
        new_theta[ncolors + k] = sigma2;

        // Posterior for mu_k: N(mu_n, sigma_n^2) with ordering
        // sigma_n^2 = 1 / (1/sigma2_0 + n_k/sigma_k^2)
        // mu_n = sigma_n^2 * (mu_0/sigma2_0 + sum_y_k/sigma_k^2)
        double prec_prior = 1.0 / sigma2_0;
        double prec_lik = n_k / sigma2;
        double sigma2_n = 1.0 / (prec_prior + prec_lik);
        double mu_n = sigma2_n * (prec_prior * mu_0 + sum_y[k] / sigma2);
        double mu_sd = std::sqrt(sigma2_n);

        double lo = (k > 0) ? new_theta[k - 1] : -1e10;
        double hi = (k + 1 < ncolors) ? new_theta[k + 1] : 1e10;
        new_theta[k] = rng.truncated_normal(mu_n, mu_sd, lo, hi);
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
        double upper = (k + 1 < ncolors) ? new_theta[k + 1] : 1e10;
        new_theta[k] = rng.truncated_gamma(alpha_n, beta_n, lower, upper);
        lower = new_theta[k];
      }
      break;
    }
    default:
      throw std::invalid_argument("Unsupported emission family type");
  }

  return new_theta;
}

}  // namespace mdgm
