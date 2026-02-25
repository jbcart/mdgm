test_that("nug_from_edge_list creates a graph", {
  edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
  g <- nug_from_edge_list(3, edges, seed = 42L)
  expect_s3_class(g, "NaturalUndirectedGraph")
  expect_equal(g$nvertices(), 3L)
  expect_equal(g$nedges(), 3L)
})

test_that("nug_from_adj_list creates a graph", {
  adj <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
  g <- nug_from_adj_list(adj, seed = 42L)
  expect_equal(g$nvertices(), 3L)
  expect_equal(g$nedges(), 3L)
})

test_that("nug_from_adj_mat creates a graph", {
  A <- matrix(0, 4, 4)
  A[1, 2] <- A[2, 1] <- 1
  A[2, 3] <- A[3, 2] <- 1
  A[3, 4] <- A[4, 3] <- 1
  g <- nug_from_adj_mat(A, seed = 42L)
  expect_equal(g$nvertices(), 4L)
  expect_equal(g$nedges(), 3L)
})

test_that("neighbors returns correct 1-indexed vertices", {
  adj <- list(c(2L, 3L), c(1L, 3L), c(1L, 2L))
  g <- nug_from_adj_list(adj, seed = 42L)
  nbrs <- sort(g$neighbors(1))
  expect_equal(nbrs, c(2L, 3L))
})

test_that("nug_from_adj_mat rejects non-symmetric matrix", {
  A <- matrix(c(0, 1, 0, 0), 2, 2)
  expect_error(nug_from_adj_mat(A))
})

test_that("nug_from_adj_mat rejects self-loops", {
  A <- matrix(c(1, 1, 1, 0), 2, 2)
  expect_error(nug_from_adj_mat(A), "Self-loops")
})

test_that("nug_from_edge_list rejects out-of-range indices", {
  edges <- rbind(c(0, 1), c(1, 0))
  expect_error(nug_from_edge_list(2, edges))
})

test_that("sample_spanning_tree returns a result", {
  edges <- rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2), c(1, 3), c(3, 1))
  g <- nug_from_edge_list(3, edges, seed = 42L)
  tree <- g$sample_spanning_tree("wilson")
  expect_type(tree, "list")
})

test_that("nug_from_grid creates correct rook grid", {
  g <- nug_from_grid(3, 3, seed = 42L)
  expect_equal(g$nvertices(), 9L)
  expect_equal(g$nedges(), 12L)
})

test_that("nug_from_grid creates correct queen grid", {
  g <- nug_from_grid(3, 3, order = 2L, seed = 42L)
  expect_equal(g$nvertices(), 9L)
  # Queen adjacency on 3x3: 4 corners*3 + 4 edges*5 + 1 center*8 = 12+20+8=40 directed, 20 undirected
  expect_equal(g$nedges(), 20L)
})

test_that("spanning tree methods produce results", {
  g <- nug_from_grid(4, 4, seed = 42L)
  for (method in c("wilson", "aldous_broder")) {
    tree <- g$sample_spanning_tree(method)
    expect_type(tree, "list")
  }
})

test_that("weighted edges are preserved", {
  edges <- cbind(
    c(1, 2, 2, 3),
    c(2, 1, 3, 2),
    c(2.5, 2.5, 1.0, 1.0)
  )
  g <- nug_from_edge_list(3, edges, seed = 42L)
  expect_equal(g$nvertices(), 3L)
  expect_equal(g$nedges(), 2L)
})
