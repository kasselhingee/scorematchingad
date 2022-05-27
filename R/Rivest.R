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

Rivest_mats2theta <- function(k, A, evidx){
  theta <- c(Bingham_Amat2theta(A), k, evidx)
  return(theta)
}
Rivest_theta2mats <- function(theta){
  A <- Bingham_theta2Amat(theta[seq.int(1, length.out = length(theta) - 2)])
  k <- theta[length(theta) - 1]
  evidx <- theta[length(theta)]
  return(list(
    A = A,
    k = k,
    evidx = evidx
  ))
}

#evidx specifies the eigenvector of A that is m. Vectors are sorted increasing size
qdRivest <- function(u, k, A, evidx){
  evidx <- floor(evidx + 0.1)
  A_es <- eigen(A)
  sizeorder <- order(abs(A_es$values))
  m <- A_es$vectors[, sizeorder == evidx]
  if (m[1] > 0){m <- -m} #I think the Cpp eigen package always has the first element as negative, but I couldn't find any documentation
  qd <- exp(k * m %*% u + t(u) %*% A %*% u)
  return(qd)
}

#derivative wrt u
lldRivest_du <- function(u, k, A, evidx){
  evidx <- floor(evidx + 0.1)
  A_es <- eigen(A)
  sizeorder <- order(abs(A_es$values))
  m <- A_es$vectors[, sizeorder == evidx, drop = FALSE]
  if (m[1] > 0){m <- -m} #I think the Cpp eigen package always has the first element as negative, but I couldn't find any documentation
  du <- k * m + 2 * A %*% u
  return(du)
}
