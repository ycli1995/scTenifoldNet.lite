
#' Counts-Per-Million (CPM) Normalization
#'
#' Performs counts-per-million (CPM) normalization on a gene expression matrix.
#' This function scales each column (cell) by its total count and multiplies by
#' a constant scaling factor, typically one million, to normalize sequencing
#' depth across cells.
#'
#' @param X A gene expression matrix with genes in rows and cells in columns.
#' @param ... Additional arguments passed to methods
#'
#' @return A normalized expression matrix of the same class and dimensions.
#'
#' @export cpmNormalization
cpmNormalization <- function(X, ...) {
  UseMethod("cpmNormalization", X)
}

#' @param scale_factor Numeric scaling factor applied after column
#' normalization. Default is \code{1e6}, corresponding to standard CPM
#' normalization.
#'
#' @details
#' CPM normalization rescales each cell by its library size:
#' \deqn{
#' x_{ij}^{\mathrm{CPM}} = \frac{x_{ij}}{\sum_i x_{ij}} \times s
#' }
#' where \eqn{x_{ij}} is the count of gene \eqn{i} in cell \eqn{j}, and
#' \eqn{s} is the scaling factor (default \eqn{10^6}).
#'
#' This normalization removes differences in sequencing depth across cells while
#' preserving relative expression abundance.
#'
#' @examples
#' # Dense matrix
#' X <- matrix(c(10, 20, 30, 40), nrow = 2)
#' colnames(X) <- c("Cell1", "Cell2")
#' rownames(X) <- c("Gene1", "Gene2")
#' cpmNormalization(X)
#'
#' @rdname cpmNormalization
#' @export
#' @method cpmNormalization matrix
cpmNormalization.matrix <- function(X, scale_factor = 1e6, ...) {
  old.dimnames <- dimnames(X)
  X <- matrix_rc_norm(X, scale_factor)
  dimnames(X) <- old.dimnames
  X
}

#' @examples
#' # Sparse matrix
#' Xs <- Matrix::Matrix(X, sparse = TRUE)
#' cpmNormalization(Xs)
#'
#' @rdname cpmNormalization
#' @export
#' @method cpmNormalization dgCMatrix
cpmNormalization.dgCMatrix <- function(X, scale_factor = 1e6, ...) {
  X@x <- X@x / rep.int(Matrix::colSums(X), diff(X@p)) * scale_factor
  X
}
