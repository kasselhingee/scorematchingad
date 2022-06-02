#' @title Rivest (1984) Model Fitting
#' @description Rivest (1984) proposed a 6-parameter version of the Fisher-Bingham distribution that
#' extended the Kent model to so that the vector `m` could have a non-zero eigenvalue in the `A` matrix.
#' @details
#' This distribution is also mentioned by Mardia and Jupp (p177, 2000).
#' There may be some instability when the first element of the `m` vector is close to zero due to the
#' method used to account for the sign-invariance of eigenvectors.
#' This method fits the matrix `A` and finds the best eigenvector for `m`.
#' To select the correct *sign* of the eigenvector `m`, the method assumes that the first element of `m` is positive,
#' and that `k` can be negative.
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

Rivest_theta2FBtheta <- function(theta){
  mats <- Rivest_theta2mats(theta)
  A_es <- eigen(mats$A)
  # print(A_es)
  evalorder <- order(abs(A_es$values), decreasing = FALSE)
  m <- A_es$vectors[, evalorder[mats$evidx]]
  # print(m)
  if (is.na(m[1])){browser()}
  if (isTRUE(m[1] < 0)){m <- -m}
  m <- sign(k) * m
  k <- abs(k)
  FB_mats2theta(k, m, mats$A)
}

#evidx specifies the eigenvector of A that is m. Vectors are sorted increasing size
qdRivest <- function(u, k, A, evidx){
  evidx <- floor(evidx + 0.1)
  A_es <- eigen(A)
  evalorder <- order(abs(A_es$values), decreasing = FALSE)
  m <- A_es$vectors[, evalorder[evidx]]
  if (m[1] < 0){m <- -m}
  # cat("m value: ", m)
  qd <- exp(k * m %*% u + t(u) %*% A %*% u)
  return(qd)
}

#derivative wrt u
lldRivest_du <- function(u, k, A, evidx){
  evidx <- floor(evidx + 0.1)
  A_es <- eigen(A)
  evalorder <- order(abs(A_es$values), decreasing = FALSE)
  m <- A_es$vectors[, evalorder[evidx], drop = FALSE]
  if (m[1] < 0){m <- -m}
  du <- k * m + 2 * A %*% u
  return(du)
}
#double derivative wrt u
lldRivest_dudu <- function(u, k, A, evidx){
  du <- 2 * A
  return(du)
}
