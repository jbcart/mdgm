#pragma once

#include <random>
#include <span>
#include <vector>

namespace mdgm {
  
class RNG {
 public:
  using rng_type = std::mt19937_64;
  using result_type = rng_type::result_type; 

  RNG() {
    reseed_random();
  }

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

  std::size_t discrete(const std::vector<double>& weights) {
    std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
    return dist(rng_);
  }

  std::size_t discrete(std::span<const double> weights) {
    std::discrete_distribution<std::size_t> dist(weights.data(), weights.data() + weights.size());
    return dist(rng_);
  }

 private:
  rng_type rng_;

  void reseed_random() {
    std::random_device rd;
    std::seed_seq seed_seq{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    rng_.seed(seed_seq);
  }
};

} // namespace mdgm
