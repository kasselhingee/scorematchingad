#' @noRd
#' @title Obtain AL and bL from Astar (the reverse is not working)
#' @description Equation 3 of the JASA paper parameterises the quadratic term of the density as a matrix A_L and a vector b_L.
#' The function `toAstar` attempted to combine these into the single matrix of parameters, Astar, as defined in equation 2.
#' @param ALs The block of A that is non zero (which is the first p-1th row and column)
#' @param bL The first p-1th elements of b (the final element of b must be zero by definition)
#' @examples
#' Astar <- rWishart(...)
#' Qin <- orthogmatwith111vec()
#' Astar <- Qin %*% diag(c(3, -2, 0)) %*% t(Qin)
#' fromAstar(Astar)
#'
#' AL = matrix(c(-166, 117, 117, -333), nrow = 2, ncol = 2)
#' bL = rep(1, 2)
#' Astar <- toAstar(AL, bL)
#' fromAstar(Astar)
#'
# @export
toAstar <- function(AL, bL){
  stop("toAstar() is faulty - it has the same result for measurements on the simplex, but not the appropriate eigenvectors.")
  p <- ncol(AL) + 1
  eigenspace <- eigen(AL)
  vL <- eigenspace$vectors
  eigenvalues <- rep(0, p)
  eigenvalues[-p] <- eigenspace$values

  Q <- matrix(NA, p, p)

  #from bL get Q[p, -p]
  Q[p,-p] <- diag(1/(2 * eigenvalues[-p])) %*% t(vL) %*% bL

  Q[-p, -p] <- vL + matrix(Q[p, -p], ncol = p-1, nrow = p - 1, byrow = TRUE)
  Q[, p] <- 1/sqrt(p)
  stopifnot(isTRUE(all.equal(sqrt(colSums(Q^2)), rep(1, p))))
  Astar <- Q %*% diag(eigenvalues) %*% t(Q)


  #check
  const <- sum(eigenvalues[-p] * Q[p, -p]^2)
  u <- runif(p)
  u <- u / sum(u)
  uL <- u[-p]
  stopifnot(abs((u %*% Astar %*% u) - (uL %*% AL %*% uL + t(bL) %*% uL + const)) < 1E-10)

  stopifnot(isSymmetric(Astar))
  vecs <- eigen(Astar)$vectors
  stopifnot(isTRUE(all.equal(vecs %*% t(vecs), diag(1, nrow = p))))
  vecs[, which.min(abs(eigen(Astar)$value))]
  return(Astar)
}

# @describeIn toAstar The reverse of toAstar()
# @export
fromAstar <- function(Astar){
  p = ncol(Astar)
  DL = Astar[1:(p-1), 1:(p-1)]
  DB = Astar[1:(p-1), p]
  DC = Astar[p, p]

  ones <- rep(1, p)
  onesL <- rep(1, p-1)

  AL <- DL - DB %*% t(onesL) - onesL %*% t(DB) + DC * onesL %*% t(onesL)
  bL <- 2 * (DB - DC * onesL)
  return(list(
    AL = AL,
    bL = bL,
    const = DC))
}
