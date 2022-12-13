# estimates of the variance of the estimator for cppad_search
sme_estvar <- function(smofun, estimate, Y, Yapproxcentres = NA * Y, approxorder = 10){
  # generate tapes with respect to the measurement (parameter is dynamic)
  p <- ncol(Y)
  Jsmofun <- pTapeJacobian(smofun, attr(smofun, "xtape"), attr(smofun, "dyntape"))
  Hsmofun <- pTapeJacobian(Jsmofun, attr(smofun, "xtape"), attr(smofun, "dyntape"))
  
  Jsmofun_u <- swapDynamic(Jsmofun, attr(smofun, "dyntape"), attr(smofun, "xtape"))
  Hsmofun_u <- swapDynamic(Hsmofun, attr(smofun, "dyntape"), attr(smofun, "xtape"))

  hess <- tape_eval_wsum(Hsmofun_u, xmat = Y, pmat = 0*estimate, xcentres = Yapproxcentres, approxorder = 10)
  sens <- -matrix(hess, byrow = TRUE, ncol = length(estimate))/nrow(Y)
  sensinv <- solve(sens)

  grads <- tape_eval(Jsmofun_u, xmat = Y, pmat = estimate, xcentres = Yapproxcentres, approxorder = 10)
  variability <- cov(grads)

  Ginfinv <- sensinv %*% variability %*% sensinv /nrow(Y) #inverse of the Godambe information matrix, also called the sandwich information matrix
  return(Ginfinv) 
}

