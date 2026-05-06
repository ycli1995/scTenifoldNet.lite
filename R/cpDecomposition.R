
#' CP Decomposition for a 3-Mode Tensor
#'
#' Performs CANDECOMP/PARAFAC (CP) decomposition on a 3-mode tensor represented
#' as a list of matrix slices. This implementation is optimized and only for the
#' 3-mode tensor case, where a tensor is represented as \eqn{K} frontal slices
#' of dimension \eqn{I \times J}.
#'
#' @param tnsr A list of numeric matrices representing frontal slices of a
#' 3-mode tensor. All matrices must have identical dimensions.
#' @param R Integer. Target CP rank (number of components).
#' @param U_list Optional list of initialized factor matrices. Must contain
#' exactly three matrices with dimensions \code{I x R}, \code{J x R}, and
#' \code{K x R}, respectively. If \code{NULL}, random initialization is used.
#' @param max_iter Integer. Maximum number of ALS iterations. Default is 25.
#' @param tol Numeric. Convergence tolerance based on relative change in
#' residual Frobenius norm. Default is \code{1e-5}.
#' @param seed Optional integer random seed for reproducible initialization.
#' Ignored if \code{U_list} is provided.
#' @param verbose Logical. Whether to display a progress bar. Default is
#' \code{TRUE}.
#'
#' @return An invisible list containing:
#' \describe{
#' \item{U}{List of 3 factor matrices corresponding to the 3 tensor modes.}
#' \item{lambdas}{Numeric vector of length \code{R} containing component
#' weights.}
#' \item{conv}{Logical indicating whether ALS converged before reaching
#' \code{max_iter}.}
#' \item{norm_percent}{Percentage of tensor Frobenius norm explained by
#' the decomposition.}
#' \item{fnorm_resid}{Final residual Frobenius norm.}
#' \item{all_resids}{Residual Frobenius norm at each ALS iteration.}
#' \item{est}{A list containing tensor metadata: \describe{
#' \item{modes}{Integer vector of tensor dimensions \code{c(I, J, K)}.}
#' \item{num_modes}{Number of tensor modes (always 3).}
#' }}
#' }
#'
#' @details
#' The function applies alternating least squares (ALS) to estimate a
#' rank-\eqn{R} CP decomposition:
#' \deqn{
#' \mathcal{X} \approx \sum_{r=1}^{R} \lambda_r ,
#' \mathbf{u}_r^{(1)} \circ \mathbf{u}_r^{(2)} \circ \mathbf{u}_r^{(3)}
#' }
#' where \eqn{\lambda_r} are component weights and \eqn{\mathbf{u}_r^{(n)}} are
#' factor vectors for mode \eqn{n}.
#'
#' To improve performance and reduce memory usage, this implementation:
#' \itemize{
#' \item avoids explicit construction of large Khatri-Rao products,
#' \item uses slice-wise matricized tensor times Khatri-Rao product (MTTKRP),
#' \item computes convergence algebraically without reconstructing the full
#' tensor.
#' }
#'
#' This function is a lightweight and optimized CP-ALS implementation. Unlike
#' general-purpose tensor decomposition libraries, it is specialized for 3-mode
#' tensors stored as lists of dense matrices, which substantially reduces memory
#' overhead in typical single-cell network applications.
#'
#' Convergence is assessed using the relative change in residual Frobenius norm:
#' \deqn{
#' \frac{||X - \hat{X}|_F^{(t)} - |X - \hat{X}|_F^{(t-1)}|}{|X|_F} < \mathrm{tol}
#' }
#'
#' @references
#' T. Kolda, B. Bader. (2009). Tensor Decompositions and Applications.
#' \emph{SIAM Review}, 51(3), 455--500.
#' doi:10.1137/07070111X
#'
#' @examples
#' set.seed(1)
#' tnsr <- replicate(4, matrix(rnorm(25), 5, 5), simplify = FALSE)
#' fit <- cpDecomposition(tnsr, R = 2, max_iter = 10, verbose = FALSE)
#'
#' fit$conv
#' fit$lambdas
#' lapply(fit$U, dim)
#'
#' @importFrom stats rnorm
#' @importFrom utils tail
#' @export
cpDecomposition <- function(
    tnsr,
    R,
    U_list = NULL,
    max_iter = 25,
    tol = 1e-5,
    seed = 42,
    verbose = TRUE
) {
  stopifnot(is.list(tnsr))
  each_dim <- unique(lapply(tnsr, dim))
  if (length(each_dim) > 1) {
    stop("The tensor slices have different dimensions.")
  }
  if (lengths(each_dim) != 2) {
    stop("The tensor slices must be 2D matrices.")
  }
  modes <- c(each_dim[[1]], length(tnsr))
  num_modes <- length(modes)
  if (length(U_list) == 0) {
    U_list <- list()
    if (length(seed) > 0) {
      set.seed(seed)
    }
    for (i in seq_along(modes)) {
      U_list[[i]] <- matrix(rnorm(modes[i] * R), nrow = modes[i], ncol = R)
    }
  }
  if (length(U_list) != num_modes) {
    stop("The random initiated `U_list` must contain ", num_modes, " matrices.")
  }
  for (i in seq_along(U_list)) {
    stopifnot(all(dim(U_list[[i]]) == c(modes[i], R)))
  }
  tnsr_norm_sq <- sum(vapply(
    X = tnsr,
    FUN = function(x) sum(x ^ 2),
    FUN.VALUE = numeric(1L)
  ))
  tnsr_norm <- sqrt(tnsr_norm_sq)

  curr_iter <- 1
  converged <- FALSE
  fnorm_resid <- rep(0, max_iter)
  lambdas <- numeric(R)
  prev_resid <- Inf
  if (verbose) {
    pb <- cli::cli_progress_bar("CP decomposition", total = max_iter)
  }
  # ================================================================
  # FAST PATH: 3-mode tensor (the only case used in scTenifoldNet)
  # Uses slicewise MTTKRP (avoids forming large Khatri-Rao products)
  # and algebraic convergence check (avoids full tensor reconstruction)
  # ================================================================
  I <- modes[1]
  J <- modes[2]
  K <- modes[3]

  slices <- tnsr

  for (curr_iter in seq_len(max_iter)) {
    if (verbose) {
      cli::cli_progress_update(id = pb, set = curr_iter)
    }
    # --- Mode 1 update ---
    V <- crossprod(U_list[[2]]) * crossprod(U_list[[3]])
    mttkrp <- matrix(0, I, R)
    for (k in seq_len(K)) {
      mttkrp <- mttkrp +
        slices[[k]] %*% U_list[[2]] * rep(U_list[[3]][k, ], each = I)
    }
    tmp <- mttkrp %*% solve(V)
    lambdas <- colSums(abs(tmp))
    U_list[[1]] <- tmp %*% diag(1/lambdas, length(lambdas))

    # --- Mode 2 update ---
    V <- crossprod(U_list[[1]]) * crossprod(U_list[[3]])
    mttkrp <- matrix(0, J, R)
    for (k in seq_len(K)) {
      mttkrp <- mttkrp +
        crossprod(slices[[k]], U_list[[1]]) * rep(U_list[[3]][k, ], each = J)
    }
    tmp <- mttkrp %*% solve(V)
    lambdas <- colSums(abs(tmp))
    U_list[[2]] <- tmp %*% diag(1/lambdas, length(lambdas))

    # --- Mode 3 update ---
    V <- crossprod(U_list[[1]]) * crossprod(U_list[[2]])
    mttkrp <- matrix(0, K, R)
    for (k in seq_len(K)) {
      mttkrp[k, ] <- colSums(U_list[[1]] * (slices[[k]] %*% U_list[[2]]))
    }
    tmp <- mttkrp %*% solve(V)
    lambdas <- colSums(abs(tmp))
    U_list[[3]] <- tmp %*% diag(1/lambdas, length(lambdas))

    # Algebraic convergence: ||X-est||^2 = ||X||^2 - 2<X,est> + ||est||^2
    inner <- sum(lambdas * colSums(U_list[[3]] * mttkrp))
    Gamma <- crossprod(U_list[[1]]) *
      crossprod(U_list[[2]]) *
      crossprod(U_list[[3]])
    est_norm_sq <- sum(outer(lambdas, lambdas) * Gamma)
    curr_resid <- sqrt(max(tnsr_norm_sq - 2 * inner + est_norm_sq, 0))

    fnorm_resid[curr_iter] <- curr_resid
    if (curr_iter > 1 && abs(curr_resid - prev_resid) / tnsr_norm < tol) {
      converged <- TRUE
      if (verbose) {
        cli::cli_progress_update(id = pb, set = max_iter, force = TRUE)
      }
      break
    }
    prev_resid <- curr_resid
  }

  if (!converged) {
    if (verbose) {
      cli::cli_progress_update(id = pb, set = max_iter, force = TRUE)
    }
  }
  if (verbose) {
    cli::cli_progress_done(id = pb)
  }

  fnorm_resid <- fnorm_resid[fnorm_resid != 0]
  norm_percent <- (1 - (tail(fnorm_resid, 1) / tnsr_norm)) * 100
  est <- list(modes = modes, num_modes = num_modes)
  invisible(list(
    U = U_list,
    lambdas = lambdas,
    conv = converged,
    norm_percent = norm_percent,
    fnorm_resid = tail(fnorm_resid, 1),
    all_resids = fnorm_resid,
    est = est
  ))
}
