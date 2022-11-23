#' @noRd
#' @title Draft standard error estimates for a CppAD-based estimator
#' @description The standard errors of an estimator via the Godambe information matrix (or sandwich information matrix).
#' @param smofun A tape of the score matching objective calculation
#' @param est The parameter set
#' @param utabl A matrix of observations, each row being an observation.
#' @param Other Other arguments are the same as in `smobj_sum()`, `smobjgrad_sum()` and `smobjhess_sum()`. 
#' @details It seems likely that `sqrt(n) * (est - true)` will have an asympotically normal distribution with zero mean and a covariance of
#' given by the inverse of the Godambe information matrix, `invG`.
#' 
#' The weights `w` are assumed fix so that the minimiser of the score-matching divergence is still an M-estimator.
#' The sensitivity matrix is estimated as
#' \deqn{\hat{H(\theta)} = \hat{E} -grad(grad(w smo(\theta;Y))),}
#' the variability matrix is estimated as
#' \deqn{\hat{J}(\theta) = var(grad(w smo(\theta;Y))),}
#' and the variance of estimator is estimated as
#' \deqn{\hat{H}(\theta)^{-1}\hat{J}(\theta)\hat{H}(\theta)^{-1}/n,}
#' where `n` is the number of observations, regardless of their weight.
#' @return The square root of the diagonal of the matrix `invG/n`, which is an estimate of estimator SE.
# @export
cppadSE <- function(smofun, theta, utabl, ...){
  sqrt(diag(cppadSEsq(smofun, theta, utabl, ...)))
}

# @describeIn cppadSE Returns the matrix `invG/n`.
cppadSEsq <- function(smofun, theta, utabl,
                      Jsmofun_u = NULL, Hsmofun_u = NULL,
                      uboundary = NULL, boundaryapprox = NULL, approxorder = NULL, w = NULL, wboundary = NULL){


  # convert old-style boundary specification into new style
  boundaryapprox <- rbind(NA * utabl, boundaryapprox)
  utabl <- rbind(utabl, uboundary)
  w <- c(w, wboundary)
  if (is.null(w)){w <- rep(1, nrow(utabl))}
  stopifnot(all(w[[1]] == w))

  theta <- as.vector(theta)

  Ginfinv <- sme_estvar(smofun, theta, Y = utabl, Yapproxcentres = boundaryapprox, approxorder = 10)
  return(Ginfinv) 
}

# estimates of the variance of the estimator
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

