#include <cpp11.hpp>
#include <mdgm/rng.hpp>
#include <memory>

[[cpp11::register]]
cpp11::external_pointer<mdgm::RNG> rng_create_cpp() {
  auto rng = std::make_unique<mdgm::RNG>();
  return cpp11::external_pointer<mdgm::RNG>(rng.release());
}

[[cpp11::register]]
cpp11::external_pointer<mdgm::RNG> rng_create_seed_cpp(int seed) {
  auto rng = std::make_unique<mdgm::RNG>(seed);
  return cpp11::external_pointer<mdgm::RNG>(rng.release());
}

[[cpp11::register]]
void rng_set_seed_cpp(cpp11::external_pointer<mdgm::RNG> rng, int seed) {
  rng->Reseed(static_cast<mdgm::RNG::result_type>(seed));
}

[[cpp11::register]]
void rng_reseed_random_cpp(cpp11::external_pointer<mdgm::RNG> rng) {
  rng->Reseed();
}

[[cpp11::register]]
double rng_uniform_cpp(cpp11::external_pointer<mdgm::RNG> rng, double a, double b) {
  return rng->uniform(a, b);
}

[[cpp11::register]]
int rng_uniform_int_cpp(cpp11::external_pointer<mdgm::RNG> rng, int a, int b) {
  return rng->uniform(a, b);
}

[[cpp11::register]]
double rng_normal_cpp(cpp11::external_pointer<mdgm::RNG> rng, double mean, double stddev) {
  return rng->normal(mean, stddev);
}

[[cpp11::register]]
int rng_discrete_cpp(cpp11::external_pointer<mdgm::RNG> rng, const cpp11::doubles& weights) {
  std::vector<double> w_vec(weights.begin(), weights.end());
  return static_cast<int>(rng->discrete(w_vec));
}

