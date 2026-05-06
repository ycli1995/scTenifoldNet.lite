
#' Identify Differentially Regulated Genes from Aligned Manifolds
#'
#' Quantifies differential gene regulation between two aligned gene regulatory
#' networks by measuring the distance between matched gene embeddings in a
#' shared manifold space.
#'
#' @param manifoldOutput A numeric matrix returned by
#' \code{\link{manifoldAlignment}}, containing aligned manifold coordinates.
#' Row names must follow the format \code{"X_gene"} and \code{"Y_gene"},
#' with genes from network \code{X} followed by the corresponding genes from
#' network \code{Y} in the same order.
#' @param gKO Optional character vector of genes to exclude when estimating the
#' background variance used for fold-change calculation. This is typically
#' used to exclude known perturbed genes (e.g., knockout targets) from the
#' null distribution. Default is \code{NULL}.
#' @param verbose Logical. Whether to print progress and summary messages.
#' Default is \code{TRUE}.
#'
#' @return A data frame with one row per gene and the following columns:
#' \describe{
#' \item{gene}{Gene name.}
#' \item{distance}{Euclidean distance between aligned \code{X} and \code{Y}
#' embeddings for the gene.}
#' \item{Z}{Z-score of the Box-Cox transformed distance.}
#' \item{FC}{Squared distance normalized by the background mean squared
#' distance.}
#' \item{p.value}{Chi-squared test p-value for differential regulation.}
#' \item{p.adj}{False discovery rate (FDR)-adjusted p-value.}
#' }
#' Rows are ordered by decreasing \code{distance}, such that the most
#' differentially regulated genes appear first.
#'
#' @details
#' This function takes the output of \code{\link{manifoldAlignment}} and
#' computes a per-gene differential regulation score based on the Euclidean
#' distance between each gene's coordinates in the two aligned manifolds:
#' \deqn{
#' d_i = |X_i - Y_i|2
#' }
#' where \eqn{X_i} and \eqn{Y_i} are the aligned coordinates of gene
#' \eqn{i} in the two manifolds.
#'
#' To stabilize variance and improve normality, distances are optionally
#' transformed using a Box-Cox power transformation. The transformed distances
#' are then standardized to Z-scores.
#'
#' A fold-change-like statistic is computed as:
#' \deqn{
#' FC_i = \frac{d_i^2}{E[d^2]}
#' }
#' where the denominator is the mean squared distance across background genes
#' (excluding \code{gKO}, if provided).
#'
#' Statistical significance is assessed using a chi-squared approximation with
#' one degree of freedom, followed by false discovery rate correction.
#'
#' Genes with large distances and small adjusted p-values are interpreted as
#' having significantly altered regulatory relationships between the two
#' conditions.
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
#'
#' emb <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
#' rownames(emb) <- c(paste0("X_", genes), paste0("Y_", genes))
#'
#' res <- dRegulation(emb, verbose = FALSE)
#'
#' head(res)
#'
#' @importFrom MASS boxcox
#' @importFrom stats dist p.adjust pchisq
#' @importFrom cli cli_alert_info cli_alert_success
#' @export
dRegulation <- function(manifoldOutput, gKO = NULL, verbose = TRUE) {
  geneList <- rownames(manifoldOutput)
  geneList <- geneList[grepl('^X_', geneList)]
  geneList <- gsub('^X_', '', geneList)
  nGenes <- length(geneList)

  eGeneList <- rownames(manifoldOutput)
  eGeneList <- eGeneList[grepl('^Y_', eGeneList)]
  eGeneList <- gsub('^Y_', '', eGeneList)
  eGenes <- length(eGeneList)

  if (nGenes != eGenes) {
    stop('Number of identified and expected genes are not the same')
  }
  if (!all(eGeneList == geneList)) {
    stop(
      'Genes are not ordered as expected. ',
      'X_ genes should be followed by Y_ genes in the same order'
    )
  }
  if (verbose) {
    cli::cli_alert_info("Computing distances for {nGenes} genes")
  }
  dMetric <- vapply(
    X = geneList,
    FUN = function(G) {
      genes <- paste0(c("X_", "Y_"), G)
      as.numeric(dist(manifoldOutput[genes, , drop = FALSE]))
    },
    FUN.VALUE = numeric(1L)
  )

  # Box-Cox transformation
  lambdaValues <- seq(-2, 2, length.out = 1000)
  lambdaValues <- lambdaValues[lambdaValues != 0]
  BC <- try(
    MASS::boxcox(dMetric ~ 1, plot = FALSE, lambda = lambdaValues),
    silent = TRUE
  )
  if (inherits(BC, 'try-error')) {
    nD <- dMetric
  } else {
    BC <- BC$x[which.max(BC$y)]
    if (BC < 0) {
      nD <- 1 / (dMetric ^ BC)
    } else {
      nD <- dMetric ^ BC
    }
  }

  Z <- scale(nD)
  E <- mean(dMetric[!geneList %in% gKO]^2)
  FC <- dMetric^2 / E
  pValues <- pchisq(q = FC, df = 1, lower.tail = FALSE)
  pAdjusted <- p.adjust(pValues, method = 'fdr')

  dOut <- data.frame(
    gene = geneList,
    distance = dMetric,
    Z = Z,
    FC = FC,
    p.value = pValues,
    p.adj = pAdjusted,
    row.names = geneList
  )
  dOut <- dOut[order(dOut$distance, decreasing = TRUE), , drop = FALSE]
  if (verbose) {
    nSig <- sum(dOut$p.adj < 0.05)
    cli::cli_alert_success(paste(
      "Differential regulation complete: ",
      "{nSig}/{nGenes} significant genes (FDR < 0.05)"
    ))
  }
  dOut
}
