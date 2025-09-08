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
invvec <- function(v, m){ matrix(v, nrow = m) }

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

# function f must have atomic inputs only
memoise_simple_function <- function(f) {
  cache <- new.env(parent = emptyenv())  # Create an isolated environment
  
  function(...) {
    key <- paste(..., collapse = "_")  # Simple key (works for atomic inputs)
    if (exists(key, envir = cache)) {
      return(get(key, envir = cache))  # Return cached result
    }
    
    result <- f(...)  # Compute function
    assign(key, result, envir = cache)  # Store result in cache
    return(result)
  }
}

commutation_mat_direct <- function(m, p){
  idxvec <- seq.int(1, m*p)
  idxmat <- invvec(idxvec, m)
  vectidxmat <- vec(t(idxmat))
  out <- diag(m*p)[match(vectidxmat, idxvec), ]
  out
}


#' @title Commutation matrix
#' @description
#' The \insertCite{magnus2019ma}{scorematchingad} commutation matrix \eqn{K_{m,p}} that is such that
#' `K %*% vec(A) == vec(t(A))` for any m x p matrix `A`.
#' @param m Number of rows of matrix `A`
#' @param p Number of columns of matrix `A`
#' @return Matrix of size mp x mp.
#' @details
#' Uses memoisation to avoid recalculating the matrix for the same set of dimensions
#' @export
commutation_mat <- memoise_simple_function(commutation_mat_direct)
