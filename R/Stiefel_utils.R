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

#' @noRd
#' @title Matrix for projection from Stiefel manifold to tangent space
#' @description
#' The projection of a matrix `Z` in the Stiefel manifold to the tangent space at `A` in the Stiefel manifold.
#' The formula is from \insertCite{@eq2.4, @edelman1998ge}{scorematchingad} and has been rearranged for vec(Z)
#' @name Stiefel_projection
NULL

#' @describeIn Stiefel_projection Function that projects `Z` onto the tangent space at `A` according to \insertCite{@eq2.4, @edelman1998ge}{scorematchingad}.
Stiefel_proj <- function(Z, A){
  (diag(nrow(A)) - 0.5 * A %*% t(A)) %*% Z - 0.5 * A %*% t(Z) %*% A
}

#' @describeIn Stiefel_projection Projection matrix for vectorised representation of `Z`. So that `invvec(Stiefel_projmat(A) %*% vec(Z)) = Stiefel_proj(Z, A)`
Stiefel_projmat <- function(A){
  diag(ncol(A)) %x% (diag(nrow(A)) - 0.5 * A %*% t(A)) - 0.5 * (t(A) %x% A) %*% commutation_mat(nrow(A), ncol(A))
}

#' @describeIn Stiefel_projection Element-wise derivative of `Stiefel_projmat()`. `i` and `j` specify the row and column index of the element respectively.
Stiefel_projmat_d <- function(A, i, j){
  zeros_m <- rep(0, nrow(A))
  zeros_p <- rep(0, ncol(A))
  ei <- zeros_m; ei[i] <- 1
  ej <- zeros_p; ej[j] <- 1
  -0.5 * (diag(ncol(A)) %x% (ei %*% t(ej) %*% t(A) + A %*% ej %*% t(ei))  + 
            ( (ej %*% t(ei)) %x% A ) %*% commutation_mat(nrow(A), ncol(A)) + 
            (t(A) %x% (ei %*% t(ej))) %*% commutation_mat(nrow(A), ncol(A)) )
}


# for use by `tape_uld_inbuilt()`
Stiefel_MF_default_xtheta <- function(amdim, x, theta){
  if (is.null(x)){stop("Default taping values for Stiefel_MF cannot be chosen using ambient dimension (amdim) alone")}
  theta <- switch(1+is.null(theta), theta, seq.int(1, amdim))
  return(list(x = x, theta = theta))
}
