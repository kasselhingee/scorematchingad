#' @title Stiefel Manifold Embedding Tools
#' @description
#' An m x p matrix `A`, p <= m, is an element of the m,p Stiefel manifold if and only `t(A) %*% A == diag(p)`.
#' A Stiefel manifold using the matrix Frobenius norm as a metric can be embedded in an ambient mp Euclidean space  \insertCite{@Section 2.2, @edelman1998ge}{scorematchingad}
#' where each matrix `A` in the Stiefel manifold is represented by a mp-length vector `vec(A)` in Euclidean space.
#' @name Stiefel_embedding
NULL

#' @describeIn Stiefel_embedding Stack columns of a matrix `A` on top of each other to obtain a vector of all elements of `A`. This wraps `as.vector()`, which does exactly the desired process.
#' @export
vec <- function(A){ as.vector(A) }

#' @describeIn Stiefel_embedding Break up a vector into columns of length `m`, and return as a matrix.
#' @export
invvec <- function(v, m){ as.matrix(v, nrow = m) }

#' @describeIn Stiefel_embedding Check is matrix is an element of Stiefel manifold.
#' @export
is_orthonormal <- function(A, tol = 1E-8){
  diff <- t(A) %*% A - diag(ncol(A))
  norm(diff, type = "F") < tol
}


#' @describeIn Stiefel_embedding Check a vector represents an element of Stiefel manifold.
#' @export
in_Stiefel <- function(v, m, tol = 1E-8){
  A <- invvec(v, m)
  is_orthonormal(A, tol = tol)
}
