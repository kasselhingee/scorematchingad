# estimates of the variance of the estimator for cppad_search
sme_estvar <- function(smofun, estimate, Y, Yapproxcentres = NA * Y, approxorder = 10){
  stopifnot(inherits(smofun, "ADFun"))
  # generate tapes with respect to the measurement (parameter is dynamic)
  p <- ncol(Y)
  Jsmofun <- tapeJacobian(smofun)
  Hsmofun <- tapeJacobian(Jsmofun)
  
  Jsmofun_u <- tapeSwap(Jsmofun)
  Hsmofun_u <- tapeSwap(Hsmofun)

  hess <- evaltape_wsum(Hsmofun_u, xmat = Y, pmat = 0*estimate, xcentres = Yapproxcentres, approxorder = 10)
  sens <- -matrix(hess, byrow = TRUE, ncol = length(estimate))/nrow(Y)
  sensinv <- solve(sens)

  grads <- evaltape(Jsmofun_u, xmat = Y, pmat = estimate, xcentres = Yapproxcentres, approxorder = 10)
  variability <- cov(grads)

  Ginfinv <- sensinv %*% variability %*% sensinv /nrow(Y) #inverse of the Godambe information matrix, also called the sandwich information matrix
  return(Ginfinv) 
}

