#' @title To or From vector form of parameters for PPI
#' @export
toPPIparamvec <- function(ALs, bL, beta0){
  theta <- c(diag(ALs), ALs[upper.tri(ALs)], bL, beta0)
  return(theta)
}

#' @describeIn toPPIparamvec Reverse of `toPPIparamvec()`.
#' @export
fromPPIparamvec <- function(theta, p){
  ALs <- matrix(NA, nrow = p - 1, ncol = p - 1)
  diag(ALs) <- theta[1:(p - 1)]
  ALs[upper.tri(ALs)] <- theta[p -1 + 1:((p-2) * (p-1)/2)]
  ALs[lower.tri(ALs)] <- ALs[upper.tri(ALs)]
  bL <- theta[p - 1 + ((p-2) * (p-1)/2) + 1:(p-1)]
  beta0 <- theta[(p - 1 + ((p-2) * (p-1)/2) + (p-1) + 1):length(theta)]
  return(list(
    ALs = ALs,
    bL = bL,
    beta0 = beta0
  ))
}
