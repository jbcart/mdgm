#' @keywords internal
#' @importFrom R6 R6Class
RNG <- R6::R6Class(
  classname = "RNG",
  cloneable = FALSE,
  public = list(
    initialize = function(seed = NULL) {
      if (!is.null(seed)) {
        private$.rng <- rng_create_seed_cpp(as.integer(seed))
      } else {
        private$.rng <- rng_create_cpp()
      }
    },

    set_seed = function(seed) {
      if (missing(seed)) {
        rng_reseed_random_cpp(private$.rng)
      } else {
        if (!is.integer(seed)) {
          seed <- as.integer(seed)
        }
        rng_set_seed_cpp(private$.rng, seed)
      }
    },

    uniform_int = function(a, b) {
      rng_uniform_int_cpp(private$.rng, as.integer(a), as.integer(b))
    },

    uniform = function(a = 0.0, b = 1.0) {
      rng_uniform_cpp(private$.rng, as.numeric(a), as.numeric(b))
    },

    discrete = function(probs) {
      rng_discrete_cpp(private$.rng, as.numeric(probs))
    },

    normal = function(mean = 0.0, sd = 1.0) {
      rng_normal_cpp(private$.rng, as.numeric(mean), as.numeric(sd))
    }
  ),
  private = list(
    .rng = NULL
  )
)
