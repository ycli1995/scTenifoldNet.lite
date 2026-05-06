#' Denoise and Aggregate Networks by Tensor Decomposition
#'
#' Performs tensor-based denoising and aggregation of multiple gene regulatory
#' networks using CANDECOMP/PARAFAC (CP) decomposition. This function takes a
#' 3-mode tensor and reconstructs a denoised consensus network from a low-rank
#' tensor approximation.
#'
#' @param xList A list of numeric square matrices representing network adjacency
#' matrices. Each matrix must have identical dimensions and matching row names
#' (gene names).
#' @param K Integer. Rank of the CP decomposition (number of latent tensor
#' components). Default is 5.
#' @param maxError Numeric. Convergence tolerance passed to
#' \code{\link{cpDecomposition}}. Default is \code{1e-5}.
#' @param maxIter Integer. Maximum number of ALS iterations passed to
#' \code{\link{cpDecomposition}}. Default is \code{1000}.
#' @param nDecimal Integer. Number of decimal places to retain in the final
#' reconstructed network. Default is 1.
#' @param sparse Logical. Whether to return the reconstructed consensus network
#' as a sparse \code{CsparseMatrix}. Default is \code{FALSE}.
#' @param seed Optional integer random seed for reproducible CP initialization.
#' Default is 42.
#' @param verbose Logical. Whether to print progress and summary messages.
#' Default is \code{TRUE}.
#'
#' @return A denoised consensus network as either:
#' \itemize{
#' \item a dense numeric matrix (default), or
#' \item a sparse \code{CsparseMatrix} if \code{sparse = TRUE}.
#' }
#'
#' The returned matrix:
#' \itemize{
#' \item has the same dimensions as the input networks,
#' \item preserves input gene names as row and column names,
#' \item is normalized to the range \code{[-1, 1]},
#' \item is rounded to \code{nDecimal} decimal places.
#' }
#'
#' @details
#' Given a list of \eqn{N} gene regulatory networks, this function stacks them
#' into a 3-mode tensor of dimension \eqn{G \times G \times N}, where \eqn{G}
#' is the number of genes. It then performs a rank-\eqn{K} CP decomposition
#' using \code{\link{cpDecomposition}} and reconstructs a denoised consensus
#' network by averaging the low-rank approximation across tensor slices.
#'
#' The resulting network captures the dominant shared regulatory structure across
#' subsampled networks while suppressing noise and unstable edges.
#'
#' In \code{scTenifoldNet.lite}, this function is used to denoise multiple
#' subsampled single-cell gene regulatory networks and combine them into a
#' stable, low-rank consensus network for downstream comparison.
#'
#' @references
#' Osorio, D., Zhong, Y., Li, G., Huang, J. Z., & Cai, J. J. (2020).
#' scTenifoldNet: A Machine Learning Workflow for Constructing and Comparing
#' Transcriptome-wide Gene Regulatory Networks from Single-Cell Data.
#' \emph{Patterns}, 1(9), 100139. doi:10.1016/j.patter.2020.100139
#'
#' @examples
#' set.seed(1)
#' genes <- paste0("G", 1:10)
#' xList <- replicate(3, {
#'   mat <- matrix(rnorm(100), 10, 10)
#'   rownames(mat) <- colnames(mat) <- genes
#'   mat
#' }, simplify = FALSE)
#'
#' net <- tensorDecomposition(xList, K = 3, maxIter = 20, verbose = FALSE)
#'
#' dim(net)
#' rownames(net)[1:5]
#'
#' @importFrom methods as
#' @export
tensorDecomposition <- function(
    xList,
    K = 5,
    maxError = 1e-5,
    maxIter = 1e3,
    nDecimal = 1,
    sparse = FALSE,
    seed = 42,
    verbose = TRUE
) {
  gene_names <- unique(unlist(lapply(xList, rownames)))
  nGenes <- length(gene_names)
  if (nGenes == 0) {
    stop('Gene names are required')
  }
  if (!all(gene_names == rownames(xList[[1]]))) {
    stop("Genes differ between each subsample.")
  }
  nNet <- length(xList)
  if (nNet < 2) {
    stop("The number of networks must > 1.")
  }
  if (verbose) {
    cli::cli_alert_info("Tensor: {nGenes} x {nGenes} x {nNet} (K={K})")
  }
  # CP decomposition for X
  tensorX <- cpDecomposition(
    tnsr = xList,
    R = K,
    max_iter = maxIter,
    tol = maxError,
    seed = seed,
    verbose = verbose
  )
  U_list <- tensorX$U
  tX <- matrix(0, nGenes, nGenes)
  for (i in seq_len(nNet)) {
    tX <- tX + tcrossprod(
      U_list[[1]] %*% diag(tensorX$lambdas * U_list[[3]][i, ]),
      U_list[[2]]
    )
  }
  tX <- tX / nNet
  tX <- tX / max(abs(tX))
  tX <- round(tX, nDecimal)
  rownames(tX) <- colnames(tX) <- gene_names
  if (verbose) {
    norm_explained <- round(tensorX$norm_percent, 1)
    cli::cli_alert_success(
      "CP decomposition complete (norm explained: {norm_explained}%)"
    )
  }
  if (sparse) {
    tX <- as(tX, 'CsparseMatrix')
  }
  tX
}
