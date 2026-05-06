test_that("cpmNormalization works", {
  # Simulating of a dataset following a negative binomial distribution with high sparcity (~67%)
  nCells = 2000
  nGenes = 100
  set.seed(1)
  X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
  X <- round(X)
  X <- matrix(X, ncol = nCells)
  rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
  X <- as.matrix(X)

  # test 1
  X1 <- cpmNormalization(X)
  expect_true(inherits(X1, "matrix"))
  expect_equal(dim(X1), c(nGenes, nCells))
  expect_true(all(round(Matrix::colSums(X1)) == 1e6))

  # Input test 2
  X <- as(X, 'CsparseMatrix')
  expect_true(inherits(X, "dgCMatrix"))
  X2 <- cpmNormalization(X)
  expect_true(inherits(X2, "dgCMatrix"))
  expect_equal(dim(X2), c(nGenes, nCells))
  expect_true(all(round(Matrix::colSums(X2)) == 1e6))
})
