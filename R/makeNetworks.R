
#' Construct Multiple Gene Regulatory Networks from Expression Matrices
#'
#' Builds a list of gene regulatory networks from multiple expression matrices
#' using principal component regression. Each input matrix is independently
#' converted into a gene regulatory network using either \code{\link{pcNet}} or
#' \code{\link{pcNetFast}}, producing one inferred network per input matrix.
#'
#' @param xList A list of gene expression matrices. Each matrix must have genes
#' in rows and samples/cells in columns, with identical gene sets and matching
#' row names across all matrices.
#' @param nComp Integer. Number of principal components used for network
#' construction. Must be at least 2 and smaller than the number of genes.
#' Default is 3.
#' @param scaleScores Logical. Whether to scale network edge weights to the
#' range \code{[-1, 1]}. Default is \code{TRUE}.
#' @param symmetric Logical. Whether to symmetrize each inferred network by
#' averaging \code{A} and \code{t(A)}. Default is \code{FALSE}.
#' @param q Numeric in \code{[0, 1)}. Quantile threshold for filtering weak
#' edges in each network. Edges below this quantile are set to zero.
#' Default is 0.95.
#' @param seed Optional integer random seed for reproducibility. Default is 42.
#' @param label Optional character string used as a prefix in progress messages.
#' Default is \code{NULL}.
#' @param fast Logical. Whether to use the fast approximation
#' \code{\link{pcNetFast}} instead of the full \code{\link{pcNet}}.
#' Default is \code{FALSE}.
#' @param nCores Optional integer specifying the number of BLAS threads to use.
#' Passed to \code{\link{pcNet}} or \code{\link{pcNetFast}}.
#' @param verbose Logical. Whether to print progress and summary messages.
#' Default is \code{TRUE}.
#'
#' @return A list of inferred gene regulatory networks, one per input matrix.
#' Each element is a square numeric matrix with genes in rows and columns,
#' representing the inferred gene regulatory network for the corresponding
#' input expression matrix.
#'
#' All returned networks:
#' \itemize{
#' \item have dimensions \code{nGenes x nGenes},
#' \item preserve gene names as row and column names,
#' \item are optionally scaled and sparsified,
#' \item are directed unless \code{symmetric = TRUE}.
#' }
#'
#' @details
#' This function applies principal component regression network inference to
#' each expression matrix in \code{xList}.
#'
#' Each matrix is transposed internally so that cells are treated as rows and
#' genes as columns, matching the expected input format of
#' \code{\link{pcNet}} and \code{\link{pcNetFast}}.
#'
#' Two network construction modes are available:
#' \itemize{
#' \item \code{fast = FALSE}: uses \code{\link{pcNet}}, the full
#' per-gene principal component regression approach.
#' \item \code{fast = TRUE}: uses \code{\link{pcNetFast}}, a faster global
#' approximation based on shared principal components.
#' }
#' This function is typically used after subsampling single-cell expression
#' matrices and before tensor-based denoising.
#'
#' @seealso \code{\link{pcNet}}, \code{\link{pcNetFast}}
#'
#' @examples
#' set.seed(1)
#' genes <- paste0("G", 1:10)
#'
#' xList <- replicate(3, {
#' mat <- matrix(rpois(200, lambda = 5), nrow = 10, ncol = 20)
#' rownames(mat) <- genes
#' mat
#' }, simplify = FALSE)
#'
#' nets <- makeNetworks(xList, nComp = 3, fast = TRUE, verbose = FALSE)
#'
#' length(nets)
#' dim(nets[[1]])
#'
#' @importFrom Matrix t
#' @importFrom cli cli_alert_info cli_alert_success cli_progress_bar
#' cli_progress_done cli_progress_update
#'
#' @export
makeNetworks <- function(
    xList,
    nComp = 3,
    scaleScores = TRUE,
    symmetric = FALSE,
    q = 0.95,
    seed = 42,
    label = NULL,
    fast = FALSE,
    nCores = NULL,
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
  if (nComp < 2 || nComp >= nGenes) {
    stop('nComp should be >= 2 and < total number of genes')
  }
  nNet <- length(xList)
  if (verbose) {
    tag <- ifelse(!is.null(label), paste0("[", label, "] "), "")
    cli::cli_alert_info("{tag}Building {nNet} gene regulatory networks")
    cli::cli_inform("PCNet - Principal Component Network Analysis")
    cli::cli_inform("Input: {ncol(xList[[1]])} samples x {nGenes} genes")
    cli::cli_inform("Parameters: nComp={nComp}, nCores={nCores}, q={q}")
    if (fast) {
      cli::cli_inform("Using `fast = TRUE`")
    }
    id <- cli::cli_progress_bar(paste0(tag, "Networks"), total = nNet)
  }
  networks <- list()
  if (fast) {
    for (i in seq_len(nNet)) {
      networks[[i]] <- pcNetFast(
        X = Matrix::t(xList[[i]]),
        nComp = nComp,
        scaleScores = scaleScores,
        symmetric = symmetric,
        q = q,
        seed = seed,
        nCores = nCores
      )
      if (verbose) {
        cli::cli_progress_update(id = id)
      }
    }
  } else {
    for (i in seq_len(nNet)) {
      networks[[i]] <- pcNet(
        X = Matrix::t(xList[[i]]),
        nComp = nComp,
        scaleScores = scaleScores,
        symmetric = symmetric,
        q = q,
        seed = seed,
        nCores = nCores
      )
      if (verbose) {
        cli::cli_progress_update(id = id)
      }
    }
  }
  if (verbose) {
    cli::cli_progress_done(id = id)
    cli::cli_alert_success("{tag}Network construction complete: {nNet} networks")
  }
  return(networks)
}
