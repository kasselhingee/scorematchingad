#' @title Simulate from a Fisher-Bingham distribution on the sphere
#' @description Applies a rotation transformation to the matrix `A` then runs `Directional::rfb()`.
#' I am not sure why this transformation is required - but estimators only work when it is applied.
#' @param n See `Directional::rfb()`
#' @param k See `Directional::rfb()`
#' @param m See `Directional::rfb()`
#' @param A See `Directional::rfb()`
#' @details
#' Simulates from the Fisher-Bingham distribution with density proportional to
#' $$exp(k <m, x> + x^T A x)$$
#' @export
rFB <- function(n, k, m, A){
  B <- Directional::rotation(c(0, 1, 0), m)
  sample <- Directional::rfb(n, k, m, -solve(B) %*% thetamats$A %*% B)
  return(sample)
}
