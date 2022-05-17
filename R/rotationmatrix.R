#' @title Create a rotation matrix in general dimensions
#' @param vec A unit vector
#' @author Janice Scealy
#' @author Kassel Hingee
#' @return Returns a matrix that rotates the unit vector in direction `vec` to the unit vector c(1, 0, ... 0)
#' @examples
#' u <- c(1, 2, 3, 4, 5, 6)
#' Rmat <- vec2northpole(u)
#' Rmat %*% u
vec2northpole <- function(vec){
  vec <- vec / sqrt(sum(vec^2))
  p <- length(vec)
  # from Janice Scealy
  H=diag(1,p)
  H[,1]=t(t(vec))
  H[1,]=t(vec)
  vec_L=t(t(vec[2:p]))
  H[2:p,2:p]=(1/(1+H[1,1]))*vec_L%*%t(vec_L)-diag(1,sum(p,-1))
  return(solve(H))
}
