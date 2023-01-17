#' @noRd
#' @title Simulate from the Rivest model
#' @description This function uses `Directional::rfb()` to simulate from the Rivest model.
#' The direction of the Fisher compenent is chosen based on the magnitude of the eigenvalues using an index.
#' @param i The eigenvector for the Fisher component is selected as the `i`th vector in increase eigenvalue magnitude.
#' E.g. `i=1` means the eigenvector with smallest eigenvalue size.
#' @param n Number of samples
#' @param k The concentration parameter for the Fisher component. Due to the sign convention for the eigenvector, this concentration parameter can be negative.
#' @param A A symmetric matrix with trace of zero.
#' @details
#' The Fisher component eigen vector is converted to have first element negative. Only the sphere is supported.
#' @examples
#' p <- 3
#' A <- rsymmetricmatrix(p, -10, 10)
#' A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint for Bingham
#' rRivest(100, -3.2, A, 2)
#' @export
rRivest <- function(n, k, A, i){
  A_es <- eigen(A)
  sizeorder <- order(abs(A_es$values))
  m <- A_es$vectors[, sizeorder == i]
  if (m[1] > 0){m <- -m}
  if (k < 0){#alter sign of m
    m <- -m
    k <- -k
  }
  sample <- rFB(n, k, m, A)
  return(sample)
}

