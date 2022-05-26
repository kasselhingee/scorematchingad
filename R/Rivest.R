#' @title Rivest (1984) Model Fitting
#' @description Rivest (1984) proposed a 6-parameter version of the Fisher-Bingham distribution that
#' extended the Kent model to so that the vector `m` could have a non-zero eigenvalue in the `A` matrix.
#' @details
#' This distribution is also mentioned by Mardia and Jupp (p177, 2000).
#' @references
#' Rivest, L.-P. (1984). On the Information Matrix for Symmetric Distributions on the Hypersphere. The Annals of Statistics, 12(3), 1085-1089. http://www.jstor.org/stable/2240982
#' Mardia, K. V., & Jupp, P. E. (2000). Directional Statistics. Wiley.
#' @examples
#' p <- 3
#' A <- rsymmetricmatrix(p, -10, 10)
#' A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint for Bingham
#' A_es <- eigen(A)
#' m <- A_es$vectors[, 1]
#' sample <- rFB(1000, 2, m, A)
#' qdRivest(sample[1, ], 2, A, 3)

#evidx specifies the eigenvector of A that is m. Vectors are sorted increasing size
qdRivest <- function(u, k, A, evidx){
  A_es <- eigen(A)
  evalorder <- order(abs(A_es$values), decreasing = FALSE)
  m <- A_es$vectors[, evalorder[evidx]]
  print(paste0("Location of evec: ", 4 - evalorder[evidx])) #the minus accounts for R's eigen() returns evals in DECREASING order while Cpp's eigen returns them in increasing
  if (m[1] > 0){m <- -m} #I think the Cpp eigen package always has the first element as negative, but I couldn't find any documentation
  print(paste0("first m value: ", m[1]))
  qd <- exp(k * m %*% u + t(u) %*% A %*% u)
  return(qd)
}

#derivative wrt u
lldRivest_du <- function(u, k, A, evidx){
  A_es <- eigen(A)
  evalorder <- order(abs(A_es$values), decreasing = FALSE)
  m <- A_es$vectors[, evalorder[evidx], drop = FALSE]
  if (m[1] > 0){m <- -m} #I think the Cpp eigen package always has the first element as negative, but I couldn't find any documentation
  du <- k * m + 2 * A %*% u
  return(du)
}
