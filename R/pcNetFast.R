
#' Construct a Gene Regulatory Network Using Fast Global Principal Component
#' Regression
#'
#' Infers a gene regulatory network from a gene expression matrix using a fast,
#' global principal component regression (PCR) approximation.
#'
#' @param X A numeric gene expression matrix with cells in rows and genes in
#' columns. Must be either a base \code{matrix} or a sparse
#' \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}. All genes must have non-zero
#' column sums.
#' @param nComp Integer. Number of global principal components used for
#' regression. Must be at least 2 and smaller than the number of genes.
#' Default is 3.
#' @param scaleScores Logical. Whether to scale the final network edge weights
#' to the range \code{[-1, 1]}. Default is \code{TRUE}.
#' @param symmetric Logical. Whether to symmetrize the inferred network by
#' averaging \code{A} and \code{t(A)}. Default is \code{FALSE}.
#' @param q Numeric in \code{[0, 1)}. Quantile threshold for filtering weak
#' edges based on absolute coefficient magnitude. Edges below this quantile
#' are set to zero. Default is 0 (no filtering).
#' @param seed Optional integer random seed for reproducible truncated SVD
#' initialization. Default is 42.
#' @param nCores Optional integer specifying the number of BLAS threads to use.
#' If \code{NULL}, half of available BLAS threads are used by default.
#'
#' @return A numeric square matrix representing the inferred gene regulatory
#' network. Rows correspond to target genes and columns correspond to
#' predictor genes. Entry \code{(i, j)} represents the inferred regulatory
#' effect of gene \code{j} on gene \code{i}.
#'
#' The returned matrix:
#' \itemize{
#' \item has dimensions \code{nGenes x nGenes},
#' \item is directed unless \code{symmetric = TRUE},
#' \item has diagonal entries forced to zero,
#' \item optionally scaled to \code{[-1, 1]},
#' \item optionally filtered to remove weak edges.
#' }
#'
#' @details
#' Unlike \code{\link{pcNet}}, which performs one regression per target gene,
#' this method computes a shared low-dimensional representation of the full
#' expression matrix and estimates all gene-gene regulatory coefficients in a
#' single matrix operation. This approximation substantially improves
#' computational efficiency and is designed for rapid network construction in
#' large single-cell transcriptomic datasets.
#'
#' \code{pcNetFast} computes a single truncated singular
#' value decomposition (SVD) of the scaled expression matrix:
#' \deqn{
#' X \approx UDV^\top
#' }
#'
#' The top \code{nComp} right singular vectors are used as global principal
#' components, and all gene-gene regression coefficients are estimated
#' simultaneously in the shared latent space. This yields a fast low-rank
#' approximation to the full PCR network:
#' \deqn{
#' W \approx X^\top Z (Z^\top Z)^{-1} V^\top
#' }
#' where \eqn{Z} is the matrix of projected principal component scores.
#'
#' Compared with \code{\link{pcNet}}, this method:
#' \itemize{
#' \item is substantially faster for large matrices,
#' \item avoids fitting one model per gene,
#' \item provides a low-rank approximation to the full PCR network,
#' \item is especially suitable for large-scale exploratory analyses.
#' }
#'
#' Because all genes share the same latent basis, \code{pcNetFast} is an
#' approximation and may be less accurate than \code{\link{pcNet}} for
#' recovering gene-specific local structure, but is typically much more
#' computationally efficient.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rpois(200, lambda = 5), nrow = 20, ncol = 10)
#' colnames(X) <- paste0("G", seq_len(ncol(X)))
#'
#' net <- pcNetFast(X, nComp = 3, scaleScores = TRUE, symmetric = FALSE)
#'
#' dim(net)
#' net[1:5, 1:5]
#' @export
pcNetFast <- function(
    X,
    nComp = 3,
    scaleScores = TRUE,
    symmetric = FALSE,
    q = 0,
    seed = 42,
    nCores = NULL
) {

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  # Check quality control: all genes must have at least one count
  if (!all(Matrix::colSums(X) > 0)) {
    stop(
      "Input matrix contains genes with zero row sums. ",
      "Please apply quality control to remove low-abundance genes."
    )
  }

  # Check input type
  input_class <- class(X)[[1]]
  valid_classes <- c("matrix", "dgCMatrix")
  if (!input_class %in% valid_classes) {
    stop("Input X must be a matrix or dgCMatrix. ", "Got: ", input_class)
  }

  # Check nComp parameter
  if (nComp < 2) {
    stop("nComp must be >= 2. Got: ", nComp)
  }
  n_genes <- ncol(X)
  if (nComp >= n_genes) {
    stop("nComp must be < number of genes (", n_genes, "). Got: ", nComp)
  }

  # ============================================================================
  # DATA PREPARATION
  # ============================================================================

  # Store gene names for later
  gene_names <- colnames(X)
  # ============================================================================
  # PRINCIPAL COMPONENT REGRESSION
  # ============================================================================

  # Set up parallel backend if using multiple cores
  nCores <- nCores %||% as.integer(RhpcBLASctl::blas_get_num_procs() / 2)
  nCores <- min(nCores, RhpcBLASctl::blas_get_num_procs())
  RhpcBLASctl::blas_set_num_threads(nCores)

  # ============================================================================
  # Get global PC scores
  # ============================================================================
  X <- scale(X)
  if (length(seed) > 0) {
    set.seed(seed)
  }
  svd_result <- RSpectra::svds(A = X, k = nComp)
  principal_components <- svd_result$v
  # Result: (n_samples x nComp) matrix of PC scores
  pc_scores <- X %*% principal_components
  # Normalize PC scores by their squared norms
  # This prevents overfitting to high-variance PCs
  score_sq_norms <- Matrix::colSums(pc_scores ^ 2)
  pc_scores_normalized <- pc_scores/ rep(score_sq_norms, each = nrow(pc_scores))

  # ============================================================================
  # Regress target gene on normalized PC scores
  # ============================================================================
  network <- tcrossprod(crossprod(X, pc_scores_normalized), principal_components)
  diag(network) <- 0

  # ============================================================================
  # POST-PROCESSING
  # ============================================================================

  # Symmetrize network if requested
  if (isTRUE(symmetric)) {
    network <- (network + Matrix::t(network)) / 2
  }

  # Scaling: normalize edge weights to [-1, 1]
  if (isTRUE(scaleScores)) {
    max_abs_value <- max(
      abs(max(network, na.rm = TRUE)),
      abs(min(network, na.rm = TRUE))
    )
    if (is.finite(max_abs_value) && max_abs_value > 0) {
      network <- network / max_abs_value
    }
  }

  # Filtering: remove weak edges below quantile threshold
  if (q > 0 && q < 1) {
    abs_network <- abs(network)
    threshold <- collapse::.quantile(abs_network, q, na.rm = TRUE)
    network[abs_network < threshold] <- 0
  }

  # Force diagonal to zero (no self-loops)
  diag(network) <- 0

  # Add gene names to rows and columns
  if (length(gene_names) > 0) {
    dimnames(network) <- list(gene_names, gene_names)
  }
  network
}
