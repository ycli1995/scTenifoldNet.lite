#' Construct and Compare Single-Cell Gene Regulatory Networks
#'
#' Runs the full \code{scTenifoldNet} workflow to construct, denoise, align,
#' and compare two single-cell gene regulatory networks from expression
#' matrices representing two biological conditions.
#'
#' This function implements the end-to-end \code{scTenifoldNet.lite} pipeline
#' for differential regulatory network analysis, including quality control,
#' normalization, network construction, tensor denoising, manifold alignment,
#' and differential regulation testing.
#'
#' @param X A gene-by-cell expression matrix for condition 1. Rows correspond
#' to genes and columns correspond to cells.
#' @param Y A gene-by-cell expression matrix for condition 2. Rows correspond
#' to genes and columns correspond to cells.
#' @param qc Logical. Whether to perform quality control filtering on both
#' expression matrices. Default is \code{TRUE}.
#' @param qc_minLibSize Integer. Minimum library size required to retain a
#' cell. Default is \code{1000}.
#' @param qc_removeOutlierCells Logical. Whether to remove low-quality outlier
#' cells during QC. Default is \code{TRUE}.
#' @param qc_minPCT Numeric. Minimum fraction of cells in which a gene must be
#' detected to be retained. Default is \code{0.05}.
#' @param qc_maxMTratio Numeric. Maximum allowed mitochondrial read fraction
#' per cell. Default is \code{0.1}.
#' @param nc_nNet Integer. Number of subsampled networks to construct per
#' condition. Default is 10.
#' @param nc_nCells Integer. Number of cells sampled per network. Default is
#' 500.
#' @param nc_nComp Integer. Number of principal components used in network
#' construction. Default is 3.
#' @param nc_symmetric Logical. Whether to symmetrize inferred networks during
#' construction. Default is \code{FALSE}.
#' @param nc_scaleScores Logical. Whether to scale network edge weights to
#' \code{[-1, 1]} during network construction. Default is \code{TRUE}.
#' @param nc_q Numeric. Quantile threshold for filtering weak edges during
#' network construction. Default is \code{0.05}.
#' @param td_K Integer. Rank used in CP tensor decomposition for denoising
#' network ensembles. Default is 3.
#' @param td_nDecimal Integer. Number of decimal places retained in the
#' reconstructed tensor networks. Default is 1.
#' @param td_maxIter Integer. Maximum number of iterations for tensor
#' decomposition. Default is \code{1000}.
#' @param td_maxError Numeric. Convergence tolerance for tensor decomposition.
#' Default is \code{1e-5}.
#' @param ma_nDim Integer. Number of dimensions used in manifold alignment.
#' Default is 30.
#' @param seed Optional integer random seed for reproducibility. Default is 42.
#' @param fast Logical. Whether to use the fast approximation
#' \code{\link{pcNetFast}} for network construction. Default is
#' \code{TRUE}.
#' @param nCores Optional integer specifying the number of BLAS/OpenMP threads
#' to use in computationally intensive steps.
#' @param verbose Logical. Whether to print progress and summary messages.
#' Default is \code{TRUE}.
#'
#' @return A named list containing:
#' \describe{
#' \item{tensorNetworks}{A list containing:
#' \describe{
#' \item{X}{Sparse denoised gene regulatory network for condition X.}
#' \item{Y}{Sparse denoised gene regulatory network for condition Y.}
#' }
#' }
#' \item{manifoldAlignment}{Aligned manifold coordinates returned by
#' \code{\link{manifoldAlignment}}.}
#' \item{diffRegulation}{A data frame of differential regulation statistics
#' returned by \code{\link{dRegulation}}.}
#' }
#'
#' @details
#' \code{scTenifoldNet} implements the complete differential network comparison
#' workflow. The pipeline consists of:
#' \enumerate{
#' \item quality control filtering of both input matrices,
#' \item CPM normalization,
#' \item gene intersection across conditions,
#' \item repeated subsampling and network construction,
#' \item tensor decomposition to denoise network ensembles,
#' \item nonlinear manifold alignment of the two denoised networks,
#' \item differential regulation analysis based on aligned manifold distance.
#' }
#'
#' The resulting differential regulation statistics quantify how strongly each
#' gene's regulatory relationships differ between the two conditions.
#'
#' This workflow is designed for scalable comparison of transcriptome-wide
#' regulatory programs from single-cell expression data.
#'
#' @references
#' Osorio, D., Zhong, Y., Li, G., Huang, J. Z., & Cai, J. J. (2020).
#' scTenifoldNet: A Machine Learning Workflow for Constructing and Comparing
#' Transcriptome-wide Gene Regulatory Networks from Single-Cell Data.
#' \emph{Patterns}, 1(9), 100139. doi:10.1016/j.patter.2020.100139
#'
#' @export
scTenifoldNet <- function(
    X, Y,
    qc = TRUE,
    qc_minLibSize = 1000,
    qc_removeOutlierCells = TRUE,
    qc_minPCT = 0.05,
    qc_maxMTratio = 0.1,
    nc_nNet = 10,
    nc_nCells = 500,
    nc_nComp = 3,
    nc_symmetric = FALSE,
    nc_scaleScores = TRUE,
    nc_q = 0.05,
    td_K = 3,
    td_nDecimal = 1,
    td_maxIter = 1e3,
    td_maxError = 1e-5,
    ma_nDim = 30,
    seed = 42,
    fast = TRUE,
    nCores = NULL,
    verbose = TRUE
) {
  if (verbose) {
    cli::cli_h1("scTenifoldNet.lite::scTenifoldNet() Pipeline")
  }
  # Step 1: Quality Control
  if (isTRUE(qc)) {
    if (verbose) {
      cli::cli_alert_info("Step 1/6: Quality control")
    }
    X <- scQC2(
      X = X,
      minLibSize = qc_minLibSize,
      removeOutlierCells = qc_removeOutlierCells,
      minPCT = qc_minPCT,
      maxMTratio = qc_maxMTratio
    )
    Y <- scQC2(
      X = Y,
      minLibSize = qc_minLibSize,
      removeOutlierCells = qc_removeOutlierCells,
      minPCT = qc_minPCT,
      maxMTratio = qc_maxMTratio
    )
  }

  # Step 2: CPM Normalization
  if (verbose) {
    cli::cli_alert_info("Step 2/6: CPM normalization")
  }
  X <- cpmNormalization(X)
  Y <- cpmNormalization(Y)

  # Step 3: Gene intersection
  sharedGenes <- intersect(rownames(X), rownames(Y))
  nGenes <- length(sharedGenes)
  if (verbose) {
    cli::cli_alert_info("Shared genes: {nGenes}")
  }
  X <- X[sharedGenes, ]
  Y <- Y[sharedGenes, ]

  # Step 4: Network construction
  if (verbose) {
    cli::cli_alert_info("Step 3/6: Building gene regulatory networks")
  }
  nc_nCellsX <- nc_nCellsY <- nc_nCells
  if (fast) {
    warning(
      "Using `fast = TRUE`. nCells for down sampling must > 75% total cells",
      immediate. = TRUE, call. = FALSE
    )
    nc_nCellsX <- max(nc_nCellsX, floor(0.75 * ncol(X)))
    nc_nCellsY <- max(nc_nCellsY, floor(0.75 * ncol(Y)))
  }
  if (verbose) {
    cli::cli_alert_info(
      "Down sample {nc_nNet} expression matrices with {nc_nCells} cells"
    )
  }
  xList <- get_down_sample_matrices(
    X = X,
    N = nc_nNet,
    nCells = nc_nCellsX,
    seed = seed
  )
  yList <- get_down_sample_matrices(
    X = Y,
    N = nc_nNet,
    nCells = nc_nCellsY,
    seed = seed
  )
  xList <- makeNetworks(
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
  yList <- makeNetworks(
    xList = yList,
    nComp = nc_nComp,
    scaleScores = nc_scaleScores,
    symmetric = nc_symmetric,
    q = nc_q,
    seed = seed,
    fast = fast,
    nCores = nCores,
    verbose = verbose
  )

  # Step 5: Tensor Decomposition
  if (verbose) {
    cli::cli_alert_info("Step 4/6: Tensor decomposition")
  }
  tX <- tensorDecomposition(
    xList = xList,
    K = td_K,
    maxError = td_maxError,
    maxIter = td_maxIter,
    nDecimal = td_nDecimal,
    seed = seed,
    verbose = verbose
  )
  tY <- tensorDecomposition(
    xList = yList,
    K = td_K,
    maxError = td_maxError,
    maxIter = td_maxIter,
    nDecimal = td_nDecimal,
    seed = seed,
    verbose = verbose
  )

  # Step 6: Manifold Alignment
  if (verbose) {
    cli::cli_alert_info("Step 5/6: Manifold alignment")
  }
  mA <- manifoldAlignment(
    X = tX,
    Y = tY,
    d = ma_nDim,
    seed = seed,
    nCores = nCores,
    verbose = verbose
  )

  # Step 7: Differential Regulation
  if (verbose) {
    cli::cli_alert_info("Step 6/6: Differential regulation analysis")
  }
  dR <- dRegulation(manifoldOutput = mA, verbose = verbose)

  # Assemble output
  outputResult <- list(
    tensorNetworks = list(
      X = as(tX, "CsparseMatrix"),
      Y = as(tY, "CsparseMatrix")
    ),
    manifoldAlignment = mA,
    diffRegulation = dR
  )
  if (verbose) {
    cli::cli_alert_success("scTenifoldNet pipeline complete")
  }
  outputResult
}
