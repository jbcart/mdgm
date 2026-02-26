test_that("sample_mrf returns valid binary coloring", {
  nug <- nug_from_grid(4, 4, seed = 42L)
  z <- sample_mrf(nug, psi = 0.5, seed = 123L)

  expect_type(z, "integer")
  expect_length(z, 16L)
  expect_true(all(z %in% c(0L, 1L)))
})

test_that("sample_mrf is reproducible with seed", {
  nug <- nug_from_grid(4, 4, seed = 42L)
  z1 <- sample_mrf(nug, psi = 0.5, seed = 99L)
  z2 <- sample_mrf(nug, psi = 0.5, seed = 99L)

  expect_identical(z1, z2)
})

test_that("sample_mrf handles multi-color (Potts)", {
  nug <- nug_from_grid(3, 3, seed = 42L)
  z <- sample_mrf(nug, psi = 0.5, n_colors = 3L, n_sweeps = 50L, seed = 42L)

  expect_type(z, "integer")
  expect_length(z, 9L)
  expect_true(all(z %in% c(0L, 1L, 2L)))
})
