#' @title Create a rotation matrix in general dimensions
#' @param vec A unit vector
#' @author Kassel Hingee
#' @description Returns a matrix that rotates the unit vector in direction `vec` to the unit vector c(1, 0, ... 0) using [`rotationmatrix()`].
#' @examples
#' u <- c(1, 2, 3, 4, 5, 6)
#' Rmat <- vec2northpole(u)
#' Rmat %*% u
vec2northpole <- function(vec){
  vec <- vec / sqrt(sum(vec^2))
  p <- length(vec)
  return(rotationmatrix(c(1, rep(0, p-1)), vec))
}

#' @title Rotation Matrix From Two Vectors
#' @description Creates a rotation matrix that rotates between two vectors. The matrix was specified by \insertCite{@Section 3.2.1, @amaral2007pi; textual}, and is such that any vector perpendicular to the two vectors is unchanged (except when the two vectors are in exactly opposite directions).
#' @param a The vector to rotate `b` to.
#' @param b A vector that will be rotate to `a`.
#' @details
#' The return matrix `Q` is a rotation such that \eqn{Qb=a}, and for any vector `z` perpendicular to both `b` and `a`, \eqn{Qz=z}.
#' In the extremely rare situation that `b` = -`a`, the \insertCite{amaral2007pi; textual} method does not apply. Instead `rotationmatrix`, rotates `b` to the south pole, applies a rotation of `pi` about the second basis vector and reverses the first rotation.

#' The same method, without the case of `b = -a` is also implemented in [`Directional::rotation()`].
#' @references \insertAllCited()
#' @examples
#' a <- c(1,2,3,4,5)
#' b <- c(0,3,4,1,2)
#' rotationmatrix(a, b) %*% b
#' @export
rotationmatrix <- function(a, b){
  stopifnot(is.vector(a))
  stopifnot(is.vector(b))
  stopifnot(length(a) == length(b))
  a <- a/sqrt(a %*% a)[[1]]
  b <- b/sqrt(b %*% b)[[1]]
  if (all(a == -b)){
    d <- length(a)
    nthpole <- c(1, rep(0, d-1))
    sthpole <- -nthpole
    halfway <- c(1/sqrt(2), 1/sqrt(2), rep(0, d-2))
    b2sthpole <- rotationmat_amaral(sthpole, b)
    sth2nth <- diag(c(-1, -1, rep(1, d-2)))
    Q <- t(b2sthpole)%*%sth2nth%*%b2sthpole
  } else {
    Q <- rotationmat_amaral(a, b)
  }
  return(Q)
}

# the exact method described in amaral et al
rotationmat_amaral  <- function(a, b){ #assumes a and b are unit vectors
  ab <- (a %*% b)[[1]]
  alpha <- acos(ab)
  c <- b - a*ab
  c <- c/sqrt(c %*% c)[[1]]
  A <- a%o%c - c%o%a
  Q = diag(length(a)) + sin(alpha)*A + (cos(alpha) - 1)*(a%o%a + c%o%c)
  return(Q)
}


