
get_down_sample_matrices <- function(X, N = 10, nCells = 500, seed = 42) {
  if (nCells > 0.9 * ncol(X)) {
    stop("`nCells` (", nCells, ") > 90% of total cells (", ncol(X), ")")
  }
  sub.seeds <- NULL
  if (length(seed) > 0) {
    set.seed(seed)
    sub.seeds <- sample(1:1e6, N)
  }
  xList <- list()
  for (i in seq_len(N)) {
    if (length(sub.seeds) > 0) {
      set.seed(sub.seeds[i])
    }
    n <- 1
    mat <- X[, sample(seq_len(ncol(X)), size = nCells), drop = FALSE]
    while (!all(Matrix::rowSums(mat) > 0)) {
      n <- n + 1
      set.seed(sub.seeds[i] * n)
      mat <- X[, sample(seq_len(ncol(X)), size = nCells), drop = FALSE]
    }
    xList[[i]] <- mat
  }
  xList
}
