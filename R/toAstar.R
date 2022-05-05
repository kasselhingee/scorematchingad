#' Convert AL and bL to Astar parametrisation
#' @description Equation 3 of the JASA paper parameterises the quadratic term of the density as a matrix A and a vector b.
#' This function combines these into the single matrix of parameters, Astar, as defined in equation 2.
#' @param ALs The block of A that is non zero (which is the first p-1th row and column)
#' @param bL The first p-1th elements of b (the final element of b must be zero by definition)
#' @examples
#' AL = matrix(c(-166, 117, 117, -333), nrow = 2, ncol = 2)
#' bL = rep(1, 2)
#' toAstar(AL, bL)
#'
#' @export
toAstar <- function(ALs, bL){
  p <- nrow(ALs) + 1
  A = matrix(0, p, p) #A as definined in eqn 3
  A[1:(p-1), 1:(p-1)] <- ALs
  Astar <- A + matrix(c(bL, 0), byrow = TRUE, nrow = p, ncol = p) # each row of A element wise addition to b, gives Astar, as defined in eqn 2
  return(Astar)
}

fromAstar <- function(Astar){
  stopifnot(isSymmetric.matrix(Astar))

}
