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
