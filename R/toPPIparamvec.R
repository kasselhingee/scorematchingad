#' @title To or From vector form of parameters for PPI
#' @export
toPPIparamvec <- function(ALs, bL, beta){
  stopifnot(isSymmetric.matrix(ALs))
  theta <- c(diag(ALs), ALs[upper.tri(ALs)], bL, beta)
  return(theta)
}

toPPIcannparam <- function(ALs, bL, beta, manifold = "sphere"){
  if (manifold == "sphere"){out <- toPPIparamvec(ALs, bL, 1 + 2 * beta)}
  else if (manifold == "simplex"){out <- toPPIparamvec(ALs, bL, beta)}
  else {stop("manifold not supported")}
  return(out)
}

#' @describeIn toPPIparamvec Reverse of `toPPIparamvec()`.
#' @export
fromPPIparamvec <- function(theta, p){
  ALs <- matrix(NA, nrow = p - 1, ncol = p - 1)
  diag(ALs) <- theta[1:(p - 1)]
  ALs[upper.tri(ALs)] <- theta[p -1 + 1:((p-2) * (p-1)/2)]
  ALs[lower.tri(ALs)] <- t(ALs)[lower.tri(ALs)]
  bL <- theta[p - 1 + ((p-2) * (p-1)/2) + 1:(p-1)]
  beta <- theta[(p - 1 + ((p-2) * (p-1)/2) + (p-1) + 1):length(theta)]
  return(list(
    ALs = ALs,
    bL = bL,
    beta = beta
  ))
}

# returns the length of the theta vector for given dimension p
ppithetalength <- function(p){
  p + #the beta
  (p-1) + #the diagonal of AL
  (p-2) * (p-1)/2 + #the upper triangle of AL
  (p-1) #the bL
}
