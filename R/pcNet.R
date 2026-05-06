
#' Construct a Gene Regulatory Network by Principal Component Regression
#'
#' Infers a gene regulatory network from a gene expression matrix using
#' principal component regression (PCR). Each gene is regressed on all other
#' genes after dimensionality reduction, and the resulting regression
#' coefficients are assembled into a directed gene-by-gene adjacency matrix.
#'
#' @param X A numeric gene expression matrix with cells in rows and genes in
#' columns. Must be either a base \code{matrix} or a sparse
#' \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}. All genes must have non-zero
#' column sums.
#' @param nComp Integer. Number of principal components used in regression.
#' Must be at least 2 and smaller than the number of genes. Default is 3.
#' @param scaleScores Logical. Whether to scale the final network edge weights
#' to the range \code{[-1, 1]}. Default is \code{TRUE}.
#' @param symmetric Logical. Whether to symmetrize the inferred network by
#' averaging \code{A} and \code{t(A)}. Default is \code{FALSE}.
#' @param q Numeric in \code{[0, 1)}. Quantile threshold for filtering weak
#' edges based on absolute coefficient magnitude. Edges below this quantile
#' are set to zero. Default is 0 (no filtering).
#' @param seed Optional integer random seed used in the internal principal
#' component computation for reproducibility. Default is 42.
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
#' For each gene, the expression profile of that gene is treated as the response
#' variable, while all remaining genes are used as predictors. Principal
#' component analysis is first applied to the predictor matrix to reduce
#' dimensionality, and regression is then performed in the reduced latent space.
#' The resulting coefficients are projected back into gene space to yield
#' inferred regulatory effects.
#'
#' This procedure is repeated for every gene, producing a full directed
#' gene-by-gene adjacency matrix.
#'
#' @references
#' Osorio, D., Zhong, Y., Li, G., Huang, J. Z., & Cai, J. J. (2020).
#' scTenifoldNet: A Machine Learning Workflow for Constructing and Comparing
#' Transcriptome-wide Gene Regulatory Networks from Single-Cell Data.
#' \emph{Patterns}, 1(9), 100139. doi:10.1016/j.patter.2020.100139
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rpois(200, lambda = 5), nrow = 20, ncol = 10)
#' colnames(X) <- paste0("G", seq_len(ncol(X)))
#'
#' net <- pcNet(X, nComp = 3, scaleScores = TRUE, symmetric = FALSE)
#'
#' dim(net)
#' net[1:5, 1:5]
#'
#' @importFrom cli cli_alert_info cli_h1 cli_inform
#' @importFrom collapse .quantile
#' @importFrom RhpcBLASctl blas_get_num_procs blas_set_num_threads
#'
#' @export
pcNet <- function(
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
  # PRINCIPAL COMPONENT REGRESSION FUNCTION
  # ============================================================================

  # ============================================================================
  # PARALLEL COMPUTATION
  # ============================================================================

  # Set up parallel backend if using multiple cores
  nCores <- nCores %||% as.integer(RhpcBLASctl::blas_get_num_procs() / 2)
  nCores <- min(nCores, RhpcBLASctl::blas_get_num_procs())
  RhpcBLASctl::blas_set_num_threads(nCores)

  # Convert list of vectors to matrix (genes x genes-1)
  # Apply regression computation to all genes
  network <- matrix(0, n_genes, n_genes)
  X <- scale(X)
  for (i in seq_len(n_genes)) {
    coeff <- compute_gene_coefficients(X = X, K = i, nComp = nComp, seed = seed)
    network[i, ] <- coeff
  }

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

# Inner function: compute regression for one gene
# Input: K = gene index (1 to n_genes)
# Output: vector of regression coefficients (length n_genes-1)
#' @importFrom RSpectra svds
#' @importFrom Matrix colSums
compute_gene_coefficients <- function(X, K, nComp = 3, seed = 42) {
  # Target gene to predict
  target_gene <- X[, K]

  # Design matrix: all genes except target (n_genes-1 features)
  design_matrix <- X[, -K]

  # Step 1: Compute truncated SVD to get principal components
  # Returns right singular vectors (loadings) for dimension reduction
  if (length(seed) > 0) {
    set.seed(seed)
  }
  svd_result <- RSpectra::svds(A = design_matrix, k = nComp)
  principal_components <- svd_result$v

  # Step 2: Project design matrix onto principal components
  # Result: (n_samples x nComp) matrix of PC scores
  pc_scores <- design_matrix %*% principal_components

  # Normalize PC scores by their squared norms
  # This prevents overfitting to high-variance PCs
  score_sq_norms <- Matrix::colSums(pc_scores ^ 2)
  pc_scores_normalized <- pc_scores/ rep(score_sq_norms, each = nrow(pc_scores))

  # Step 3: Regress target gene on normalized PC scores
  # This gives coefficients in the PC space
  pc_coefficients <- Matrix::colSums(target_gene * pc_scores_normalized)

  # Transform PC coefficients back to original gene space
  # Result: n_genes regression coefficients
  gene_coefficients <- double(ncol(X))
  gene_coefficients[-K] <- as.vector(principal_components %*% pc_coefficients)
  gene_coefficients
}
