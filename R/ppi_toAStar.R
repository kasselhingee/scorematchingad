#' @name ppi_toAstar
#' @family PPI model tools
#' @title Obtain AL and bL from Astar and Reverse
#' @description \insertCite{scealy2022sc;textual}{scorecompdir} give two parameterisations of the PPI model. The function [`ppi()`] and most other PPI functions in this package use the parameterisation with matrix \eqn{A_L} and vector \eqn{b_L}. An alternative parameterisation uses a single matrix `Astar` instead of \eqn{A_L} and \eqn{b_L}.
#' @details
#' The `Astar` parametrisation rewrites the PPI density as proportional to 
#' \deqn{\exp(u^TA^*u)\prod_{i=1}^p u_i^{\beta_i},}
#' where \eqn{A^*} (`Astar`) is a \eqn{p} by \eqn{p} matrix.
#' Because \eqn{u} lies in the simplex (in particular \eqn{\sum u_i = 1}), the density is the same regardless of the value of \eqn{1^T A^* 1}=`sum(Astar)`, where \eqn{1} is the vector of ones. Thus \eqn{A_L} and \eqn{b_L} specify \eqn{A^*} up to an additive factor. In the conversion `toAstar()`, \eqn{A^*} is returned such that \eqn{1^T A^* 1 = 0}.
#' @examples
#'  Astar <- rWishart(1, 6, diag(3))[,,1]
#'  ppi_fromAstar(Astar)
#'  ppi_toAstar(ppi_fromAstar(Astar)$AL, ppi_fromAstar(Astar)$bL)
NULL

#' @param AL The `p-1` by `p-1` \eqn{A_L} matrix
#' @param bL The vector \eqn{b_L} (length `p-1`).
#' @describeIn ppi_toAstar Takes matrix `AL` and vector `bL` and returns the equivalent `Astar` with \eqn{1^T A^* 1 = }`sum(Astar) = 0`.
#' @export
ppi_toAstar <- function(AL, bL){
  # assumes DC = 0 initially
  p <- nrow(AL) + 1
  onesL <- rep(1, p-1)
  DB <- bL / 2
  DL <- AL + DB %*% t(onesL) + onesL %*% t(DB)
  Astar <- matrix(NA, p, p)
  Astar[1:(p-1), 1:(p-1)] <- DL
  Astar[1:(p-1), p] <- DB
  Astar[p,1:(p-1)] <- t(DB)
  Astar[p,p] <- 0
  # set everything to sum to 0 i.e 1Astar1 = 0
  Astar <- Astar - sum(Astar)/(p*p)
  return(Astar)
}

#' @describeIn ppi_toAstar Takes matrix `Astar` and returns the correponding `AL` and `bL` that give the same PPI density. It also returns the constant remainder that can be ignored due to normalisation of the probability density, but is useful for checking.
#' @param Astar The \eqn{A^*} matrix (a p by p symmetric matrix)
#' @export
ppi_fromAstar <- function(Astar){
  stopifnot(isSymmetric(Astar))
  p = ncol(Astar)
  DL = Astar[1:(p-1), 1:(p-1)]
  DB = Astar[1:(p-1), p]
  DC = Astar[p, p]

  onesL <- rep(1, p-1)

  AL <- DL - DB %*% t(onesL) - onesL %*% t(DB) + DC * onesL %*% t(onesL)
  bL <- 2 * (DB - DC * onesL)
  return(list(
    AL = AL,
    bL = bL,
    const = DC))
}
