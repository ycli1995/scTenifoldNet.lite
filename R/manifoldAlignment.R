
#' Align Two Gene Regulatory Networks by Nonlinear Manifold Alignment
#'
#' Performs nonlinear manifold alignment of two gene regulatory networks to
#' project them into a shared low-dimensional latent space. This function is
#' used to align two denoised gene regulatory networks and capture their
#' topological differences in a common manifold representation.
#'
#' @param X,Y A square numeric matrix representing the first gene regulatory
#' network. Rows and columns must correspond to genes, with matching row and
#' column names.
#' @param d Integer. Number of aligned manifold dimensions to return. Default is
#' 30.
#' @param seed Optional integer random seed for reproducible eigendecomposition.
#' Default is 42.
#' @param nCores Optional integer specifying the number of BLAS/OpenMP threads
#' to use. If \code{NULL}, half of available BLAS threads are used by default.
#' @param verbose Logical. Whether to print progress and summary messages.
#' Default is \code{TRUE}.
#'
#' @return A numeric matrix of aligned manifold coordinates with
#' \code{2 * nSharedGenes} rows and \code{d} columns, where:
#' \itemize{
#' \item the first \code{nSharedGenes} rows correspond to genes from
#' network \code{X},
#' \item the next \code{nSharedGenes} rows correspond to genes from
#' network \code{Y},
#' \item columns represent aligned manifold dimensions.
#' }
#' Row names are prefixed with \code{"X_"} and \code{"Y_"} to distinguish the
#' two networks, and column names are labeled \code{"NLMA 1"}, \code{"NLMA 2"},
#' etc.
#'
#' @details
#' This function implements nonlinear manifold alignment (NLMA) for comparing
#' two gene regulatory networks.
#'
#' The two input networks are first restricted to their shared genes and
#' symmetrized. A joint graph Laplacian operator is then constructed to couple
#' the two networks through corresponding genes, encouraging matched genes to
#' remain close in the aligned latent space while preserving each network's
#' internal topology.
#'
#' Rather than explicitly constructing the full \eqn{2n \times 2n} alignment
#' matrix, this implementation defines the alignment operator implicitly and
#' applies it through matrix-vector products. This substantially reduces memory
#' usage and improves performance for large networks.
#'
#' The aligned embedding is obtained by computing the eigenvectors associated
#' with the smallest nontrivial eigenvalues of the joint alignment operator
#' using \code{RSpectra::eigs()}.
#'
#' Genes with large differences between their aligned coordinates in the two
#' manifolds are interpreted as having altered regulatory relationships between
#' the two conditions.
#'
#' @references
#' Osorio, D., Zhong, Y., Li, G., Huang, J. Z., & Cai, J. J. (2020).
#' scTenifoldNet: A Machine Learning Workflow for Constructing and Comparing
#' Transcriptome-wide Gene Regulatory Networks from Single-Cell Data.
#' \emph{Patterns}, 1(9), 100139. doi:10.1016/j.patter.2020.100139
#'
#' @examples
#' set.seed(1)
#' genes <- paste0("G", 1:20)
#'
#' X <- matrix(rnorm(400), 20, 20)
#' Y <- matrix(rnorm(400), 20, 20)
#' rownames(X) <- colnames(X) <- genes
#' rownames(Y) <- colnames(Y) <- genes
#'
#' X <- (X + t(X)) / 2
#' Y <- (Y + t(Y)) / 2
#'
#' emb <- manifoldAlignment(X, Y, d = 5, verbose = FALSE)
#'
#' dim(emb)
#' rownames(emb)[1:4]
#'
#' @importFrom RhpcBLASctl blas_get_num_procs blas_set_num_threads
#' omp_set_num_threads
#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom RSpectra eigs
#' @importFrom Matrix rowSums
#'
#' @export
manifoldAlignment <- function(
    X, Y,
    d = 30,
    seed = 42,
    nCores = NULL,
    verbose = TRUE
) {
  sharedGenes <- intersect(rownames(X), rownames(Y))
  n <- length(sharedGenes)
  if (verbose) {
    cli::cli_alert_info("Manifold alignment: {n} shared genes, d={d}")
  }
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]

  X <- (X + Matrix::t(X)) / 2
  Y <- (Y + Matrix::t(Y)) / 2

  # The old construction of W is too expensive.
  # wX <- X + 1
  # wY <- Y + 1
  # mu <- 0.9 * (sum(wX) + sum(wY)) / (2 * n)
  #
  # W <- matrix(0, 2 * n, 2 * n)
  # W[1:n, 1:n] <- wX
  # W[(n + 1):(2 * n), (n + 1):(2 * n)] <- wY
  # offIdx <- seq_len(n)
  # W[cbind(offIdx, offIdx + n)] <- mu
  # W[cbind(offIdx + n, offIdx)] <- mu
  # diag(W) <- 0
  # diag(W) <- colSums(W)
  # W <- -W
  # diag(W) <- -diag(W)

  mu <- 0.9 * (sum(X) + sum(Y) + 2 * n ^ 2) / (2 * n)
  dx <- Matrix::rowSums(X) + ncol(X) + mu
  dy <- Matrix::rowSums(Y) + ncol(Y) + mu

  op <- function(z, args) {
    # Need to return (W %*% z)
    n <- args$n
    X <- args$X
    Y <- args$Y
    dx <- args$dx
    dy <- args$dy
    mu <- args$mu

    zx <- z[1:n]
    zy <- z[(n + 1):(2 * n)]

    outx <- dx * zx - as.vector(X %*% zx) - rep.int(sum(zx), n) - mu * zy
    outy <- dy * zy - as.vector(Y %*% zy) - rep.int(sum(zy), n) - mu * zx

    return(c(outx, outy))
  }

  nCores <- nCores %||% as.integer(RhpcBLASctl::blas_get_num_procs() / 2)
  nCores <- min(nCores, RhpcBLASctl::blas_get_num_procs())
  RhpcBLASctl::omp_set_num_threads(nCores)
  RhpcBLASctl::blas_set_num_threads(nCores)

  if (length(seed) > 0) {
    set.seed(seed)
  }
  E <- suppressWarnings(RSpectra::eigs(
    A = op,
    k = d * 2,
    n = 2 * n,
    which = 'SR',
    args = list(n = n, X = X, Y = Y, dx = dx, dy = dy, mu = mu)
  ))
  E$values <- Re(E$values)
  E$vectors <- Re(E$vectors)
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[, newOrder]
  E$vectors <- E$vectors[, E$values > 1e-8]
  alignedNet <- E$vectors[, seq_len(d)]
  colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
  rownames(alignedNet) <- c(
    paste0('X_', sharedGenes),
    paste0('Y_', sharedGenes)
  )
  if (verbose) {
    cli::cli_alert_success("Manifold alignment complete: {d} dimensions")
  }
  alignedNet
}
