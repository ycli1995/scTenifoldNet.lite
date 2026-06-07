#' Simulate Gene Overexpression and Identify Differential Regulation
#'
#' Runs the \code{scTenifoldOE} workflow to simulate virtual gene overexpression
#' (OE) from single-cell expression data and identify genes with altered
#' regulatory relationships after perturbation.
#'
#' @param countMatrix A gene-by-cell count matrix containing raw or normalized
#' single-cell expression values. Rows correspond to genes and columns
#' correspond to cells.
#' @param gOE Character vector of one or more gene names to simulate as
#' over expression targets. These genes must be present in
#' \code{rownames(countMatrix)}.
#' @param fold_change A numeric scalar to multiply to the `gOE` coefficences.
#' Default is 10.
#' @param cpm Logical. Whether to apply CPM normalization prior to downstream
#' analysis. Default is \code{FALSE}.
#' @param qc Logical. Whether to perform quality control filtering on the input
#' count matrix. Default is \code{TRUE}.
#' @param qc_mtThreshold Numeric. Maximum allowed mitochondrial fraction per
#' cell for quality control filtering. Default is \code{0.1}.
#' @param qc_minLSize Integer. Minimum library size per cell for quality
#' control filtering. Default is \code{1000}.
#' @param qc_minCells Integer. Minimum number of expressing cells required to
#' retain a gene. Default is \code{25}.
#' @param nc_lambda Numeric. Directionality regularization parameter. Default
#' is 0.
#' @param nc_nNet Integer. Number of subsampled networks to construct for
#' network denoising. Default is 10.
#' @param nc_nCells Integer. Number of cells sampled per network. Default is
#' 500.
#' @param nc_nComp Integer. Number of principal components used for network
#' construction. Default is 3.
#' @param nc_scaleScores Logical. Whether to scale network edge weights to
#' \code{[-1, 1]} during network construction. Default is \code{TRUE}.
#' @param nc_symmetric Logical. Whether to symmetrize inferred networks during
#' construction. Default is \code{FALSE}.
#' @param nc_q Numeric. Quantile threshold for sparsifying weak edges during
#' network construction. Default is \code{0.9}.
#' @param td_K Integer. Rank used in tensor decomposition for denoising the
#' network ensemble. Default is 3.
#' @param td_maxIter Integer. Maximum number of iterations for tensor
#' decomposition. Default is 1000.
#' @param td_maxError Numeric. Convergence tolerance for tensor decomposition.
#' Default is \code{1e-5}.
#' @param td_nDecimal Integer. Number of decimal places retained in the
#' reconstructed tensor network. Default is 3.
#' @param ma_nDim Integer. Number of dimensions used in manifold alignment.
#' Default is 2.
#' @param fast Logical. Whether to use the fast approximation
#' \code{\link{pcNetFast}} for network construction. Default is
#' \code{TRUE}.
#' @param seed Optional integer random seed for reproducibility. Default is 42.
#' @param nCores Optional integer specifying the number of BLAS/OpenMP threads
#' to use in computationally intensive steps.
#' @param verbose Logical. Whether to print progress and summary messages.
#' Default is \code{TRUE}.
#'
#' @return A named list containing:
#' \describe{
#' \item{tensorNetworks}{A list containing:
#' \describe{
#' \item{WT}{Sparse wild-type gene regulatory network.}
#' \item{OE}{Sparse overexpression gene regulatory network.}
#' }
#' }
#' \item{manifoldAlignment}{Aligned manifold coordinates returned by
#' \code{\link{manifoldAlignment}}.}
#' \item{diffRegulation}{A data frame of differential regulation statistics
#' returned by \code{\link{dRegulation}}.}
#' }
#'
#' @details
#' This function constructs a wild-type (WT) gene regulatory network from a
#' single-cell count matrix, simulates over expression of one or more target
#' genes by multiply their outgoing regulatory effects with a fold change,
#' aligns WT and OE networks in a shared manifold space, and quantifies
#' differential regulation between the two network states.
#'
#' The procedure consists of:
#' \enumerate{
#' \item optional quality control filtering of the count matrix,
#' \item repeated cell subsampling to construct an ensemble of WT networks,
#' \item principal component regression network inference,
#' \item tensor decomposition to denoise and aggregate the network ensemble,
#' \item virtual OE simulation by removing outgoing edges from \code{gOE},
#' \item nonlinear manifold alignment of WT and OE networks,
#' \item differential regulation analysis based on aligned manifold distance.
#' }
#' Genes whose manifold positions shift substantially after overexpression are
#' interpreted as downstream targets or indirectly perturbed regulators of the
#' overexpression gene.
#'
#' @references
#' Osorio, Daniel, et al. "scTenifoldKnk: An efficient virtual overexpression tool
#' for gene function predictions via single-cell gene regulatory network
#' perturbation." \emph{Patterns} 3.3 (2022).
#'
#' @importFrom cli cli_alert_info cli_alert_success cli_h1
#'
#' @export
scTenifoldOE <- function(
    countMatrix,
    gOE,
    fold_change = 10,
    cpm = FALSE,
    qc = TRUE,
    qc_mtThreshold = 0.1,
    qc_minLSize = 1000,
    qc_minCells = 25,
    nc_lambda = 0,
    nc_nNet = 10,
    nc_nCells = 500,
    nc_nComp = 3,
    nc_scaleScores = TRUE,
    nc_symmetric = FALSE,
    nc_q = 0.9,
    td_K = 3,
    td_maxIter = 1000,
    td_maxError = 1e-05,
    td_nDecimal = 3,
    ma_nDim = 2,
    fast = TRUE,
    seed = 42,
    nCores = NULL,
    verbose = TRUE
) {
  if (verbose) {
    # Start a CLI process to report progress to the user
    gOE2 <- paste(gOE, collapse = ", ")
    cli::cli_h1("scTenifoldNet.lite::scTenifoldOE() pipeline")
    cli::cli_alert_info("Simulating genes overexpression: {gOE2}")
  }
  # Check that the requested gene to knock out is present in the input matrix
  if (any(!gOE %in% rownames(countMatrix))) {
    gOE <- paste(setdiff(gOE, rownames(countMatrix)), collapse = ", ")
    stop("The following `gOE` not found in `countMatrix`: ", gOE)
  }
  # Optional quality control: filter cells and genes
  countMatrix <- scQC(
    X = countMatrix,
    mtThreshold = qc_mtThreshold,
    minLSize = qc_minLSize,
    minCells = qc_minCells,
    qc = qc
  )
  if (verbose) {
    cli::cli_alert_success(paste(
      "Count matrix quality control applied: ",
      "retained {nrow(countMatrix)} genes and {ncol(countMatrix)} cells"
    ))
  }
  if (any(!gOE %in% rownames(countMatrix))) {
    gOE <- paste(setdiff(gOE, rownames(countMatrix)), collapse = ", ")
    stop("The following `gOE` not found in `countMatrix`: ", gOE)
  }
  .check_gKO(countMatrix, gOE)

  # Get subsample of X
  if (verbose) {
    cli::cli_alert_info(
      "Down sample {nc_nNet} expression matrices with {nc_nCells} cells"
    )
  }
  xList <- get_down_sample_matrices(
    X = countMatrix,
    N = nc_nNet,
    nCells = nc_nCells,
    seed = seed
  )
  newGenes <- rownames(xList[[1]])
  if (any(!gOE %in% newGenes)) {
    gOE <- paste(setdiff(gOE, newGenes), collapse = ", ")
    stop("`gOE` were removed due to low expression: ", gOE)
  }

  # Build an ensemble of gene regulatory networks (subsample cells, use PCR)
  WT <- makeNetworks(
    xList = xList,
    nComp = nc_nComp,
    scaleScores = nc_scaleScores,
    symmetric = nc_symmetric,
    q = nc_q,
    seed = seed,
    fast = fast,
    nCores = nCores,
    verbose = verbose
  )
  if (verbose) {
    cli::cli_alert_success(paste(
      "Network construction complete",
      "(nNet = {nc_nNet}, nCells per net = {nc_nCells})"
    ))
  }

  # Tensor decomposition (CP) to denoise / approximate the ensemble of networks
  WT <- tensorDecomposition(
    xList = WT,
    K = td_K,
    maxError = td_maxError,
    maxIter = td_maxIter,
    nDecimal = td_nDecimal,
    seed = seed,
    verbose = verbose
  )
  if (verbose) {
    cli::cli_alert_success("Tensor decomposition completed (K = {td_K})")
  }

  # Extract reconstructed network and enforce directionality
  WT <- strictDirection(WT, lambda = nc_lambda)

  # Remove self-loops
  diag(WT) <- 0

  # Transpose to have genes as rows for downstream steps
  WT <- t(WT)
  if (verbose) {
    cli::cli_alert_success("Prepared WT adjacency matrix for OE simulation")
  }

  # Simulate overexpression by zeroing outgoing edges from the OE gene
  OE <- WT
  OE[gOE, ] <- OE[gOE, ] * fold_change
  if (all(OE[gOE, ] == WT[gOE, ])) {
    stop("The WT[gOE, ] is already 0. Cannot simulate overexpression for ", gOE2)
  }
  if (verbose) {
    cli::cli_alert_success("Simulated overexpression for {gOE2}")
  }

  # Align WT and OE networks into a shared low-dimensional manifold space
  MA <- manifoldAlignment(
    X = WT,
    Y = OE,
    d = ma_nDim,
    seed = seed,
    nCores = nCores,
    verbose = verbose
  )
  if (verbose) {
    cli::cli_alert_success("Manifold alignment completed (d = {ma_nDim})")
  }

  DR <- dRegulation(manifoldOutput = MA, gKO = gOE, verbose = verbose)
  if (verbose) {
    cli::cli_alert_success("Differential regulation computed for {gOE2}")
  }

  # Prepare and return results
  outputList <- list()
  outputList$tensorNetworks$WT <- as(WT, "CsparseMatrix")
  outputList$tensorNetworks$OE <- as(OE, "CsparseMatrix")
  outputList$manifoldAlignment <- MA
  outputList$diffRegulation <- DR
  # Finish CLI process and return results
  if (verbose) {
    cli::cli_alert_success("Finished scTenifoldKnk for {gOE2}")
  }
  return(outputList)
}
