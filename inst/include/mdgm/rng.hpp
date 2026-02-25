#pragma once

#include <cmath>
#include <random>
#include <span>
#include <vector>

namespace mdgm {

class RNG {
 public:
  using rng_type = std::mt19937_64;
  using result_type = rng_type::result_type;

  RNG() { reseed_random(); }

  explicit RNG(result_type seed) : rng_(seed) {}

  template <typename T>
  T uniform(T a, T b) {
    std::uniform_int_distribution<T> dist(a, b);
    return dist(rng_);
  }

  template <typename T, typename U, typename V>
  T uniform(U a, V b) {
    return uniform<T>(static_cast<T>(a), static_cast<T>(b));
  }

  double uniform(double a = 0.0, double b = 1.0) {
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng_);
  }

  double normal(double mean = 0.0, double stddev = 1.0) {
    std::normal_distribution<double> dist(mean, stddev);
    return dist(rng_);
  }

  // Gamma(shape, scale=1) via Marsaglia and Tsang's method
  double gamma(double shape) {
    if (shape < 1.0) {
      // Gamma(a) = Gamma(a+1) * U^(1/a)
      return gamma(shape + 1.0) * std::pow(uniform(), 1.0 / shape);
    }
    double d = shape - 1.0 / 3.0;
    double c = 1.0 / std::sqrt(9.0 * d);
    while (true) {
      double x = normal();
      double v = 1.0 + c * x;
      if (v <= 0.0) continue;
      v = v * v * v;
      double u = uniform();
      if (u < 1.0 - 0.0331 * (x * x) * (x * x)) return d * v;
      if (std::log(u) < 0.5 * x * x + d * (1.0 - v + std::log(v))) return d * v;
    }
  }

  // Beta(a, b) via ratio of Gamma variates
  double beta(double a, double b) {
    double x = gamma(a);
    double y = gamma(b);
    return x / (x + y);
  }

  // Truncated Beta(a, b) on [lo, hi] via rejection sampling
  double truncated_beta(double a, double b, double lo, double hi,
                        int max_attempts = 10000) {
    for (int i = 0; i < max_attempts; ++i) {
      double sample = beta(a, b);
      if (sample >= lo && sample <= hi) return sample;
    }
    return (lo + hi) / 2.0;  // fallback
  }

  std::size_t discrete(const std::vector<double>& weights) {
    std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
    return dist(rng_);
  }

  std::size_t discrete(std::span<const double> weights) {
    std::discrete_distribution<std::size_t> dist(weights.data(),
                                                  weights.data() + weights.size());
    return dist(rng_);
  }

  std::vector<std::size_t> permutation(std::size_t n) {
    std::vector<std::size_t> perm(n);
    for (std::size_t i = 0; i < n; ++i) {
      perm[i] = i;
    }
    for (std::size_t i = n - 1; i > 0; --i) {
      std::size_t j = uniform<std::size_t>(0, i);
      std::swap(perm[i], perm[j]);
    }
    return perm;
  }

  void Reseed(result_type seed) { rng_.seed(seed); }

  void Reseed() { reseed_random(); }

 private:
  rng_type rng_;

  void reseed_random() {
    std::random_device rd;
    std::seed_seq seed_seq{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    rng_.seed(seed_seq);
  }
};

}  // namespace mdgm
