#' @title Quickly generate a symmetric matrix for examples
#' @description A simple function for generating a symmetric matrix for use in examples.
#' Distributional aspects of this simulation have not been considered.
#' Elements of the matrix are simulated from a uniform distribution.
#' @param p The desired dimension of the matrix
#' @param min The minimum of the uniform distribution.
#' @param max The maximum of the uniform distribution
#' @return A `p`x`p` symmetric matrix.
#' @examples
#' rsymmetricmatrix(5)
#' @export
rsymmetricmatrix <- function(p, min = 0, max = 1){
  A <- matrix(NA, ncol = p, nrow = p)
  A[upper.tri(A)] <- runif(sum(upper.tri(A)), max = max, min = min)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- runif(p, min = min, max = max)
  return(A)
}
