#' Convert AL and bL to Astar parametrisation
#' @description Equation 3 of the JASA paper parameterises the quadratic term of the density as a matrix A and a vector b.
#' This function combines these into the single matrix of parameters, Astar, as defined in equation 2.
#' @param ALs The block of A that is non zero (which is the first p-1th row and column)
#' @param bL The first p-1th elements of b (the final element of b must be zero by definition)
#' @examples
#' Qin <- orthogmatwith111vec()
#' Astar <- Qin %*% diag(c(3, 2, 0)) %*% t(Qin)
#' fromAstar(Astar)
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
  p = ncol(Astar)
  stopifnot(isSymmetric.matrix(Astar))
  eigenspace <- eigen(Astar, symmetric = TRUE)
  # rearrange so 0 at the end
  neworder <- order(abs(eigenspace$values), decreasing = TRUE)
  eigenvalues <- eigenspace$values[neworder]
  stopifnot( abs(eigenvalues[p]) < 1E-10 )
  Q <- eigenspace$vectors[, neworder] #Q is such that Q %*% diag(eigenvalues) %*% t(Q) = Astar
  if(!all(abs(Q[1:(p-1),p] - Q[p, p]) < 1E-10)){stop("Final eigenvector of Astar is not proportional to the vector 1,1,...,1")}

  vL <- Q[-p, -p] - matrix(Q[p, -p], ncol = p-1, nrow = p - 1, byrow = TRUE) #first (p-1) eigen vectors minus their pth element
  bL <- rowSums(2 * vL %*% diag(eigenspace$values[-p] * Q[p,-p]))
    # 2 * (eigenspace$values[1] * Q[p,1] * v1n + eigenspace$values[2] * Q[p,2] * v2n)  %*% u[1:2]
  AL <- vL %*% diag(eigenspace$values[-p]) %*% t(vL)
  const <- sum(eigenspace$values[-p] * Q[p, -p]^2)

  # check
  u <- runif(p)
  u <- u / sum(u)
  uL <- u[-p]
  stopifnot(abs((u %*% Astar %*% u) - (uL %*% AL %*% uL + const + t(bL) %*% uL)) < 1E-10)

  return(list(
    AL = AL,
    bL = bL,
    const = const))
}
