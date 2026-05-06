
strictDirection <- function(X, lambda = 1) {
  if (lambda == 0) {
    return(X)
  }
  S <- as.matrix(X)
  S[abs(S) < abs(t(S))] <- 0
  if (lambda == 1) {
    return(S)
  }
  O <- (((1 - lambda) * X) + (lambda * S))
  O <- Matrix::Matrix(O)
  return(O)
}
