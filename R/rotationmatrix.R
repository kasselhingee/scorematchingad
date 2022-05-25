#' @title INTERNAL USE ONLY Create a rotation matrix in general dimensions
#' @param vec A unit vector
#' @description The rotation operation is not unique in 3+ dimensions,
#' for consistency this function is a wrapper of `Directional::rotation()`.
#' @author Kassel Hingee
#' @return Returns a matrix that rotates the unit vector in direction `vec` to the unit vector c(1, 0, ... 0)
#' @examples
#' u <- c(1, 2, 3, 4, 5, 6)
#' Rmat <- vec2northpole(u)
#' Rmat %*% u
vec2northpole <- function(vec){
  vec <- vec / sqrt(sum(vec^2))
  northpole <- c(1, rep(0, length(vec) - 1))
  out <- Directional::rotation(vec, northpole)
  return(out)
}
