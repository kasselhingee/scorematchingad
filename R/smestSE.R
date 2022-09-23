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
#' \deqn{\hat{H(\theta) = \hat{E} -grad(grad(w smo(\theta;Y))),}
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
                      uboundary = NULL, boundaryapprox = NULL, approxorder = NULL, w = NULL){
  sens <- -smobjhess(smofun, theta, utabl,
                     Hsmofun_u = Hsmofun_u,
                     uboundary = uboundary, boundaryapprox = boundaryapprox,
                     approxorder = approxorder, w = w)
  gradsmoperpt <- smobjgrad_perpt(smofun, theta, utabl,
                  Jsmofun_u = Jsmofun_u,
                  uboundary = uboundary, boundaryapprox = boundaryapprox,
                  approxorder = approxorder)
  gradsmoperpt <- do.call(rbind, gradsmoperpt)
  if (!is.null(w)) {
     gradsmoperpt <- w * gradsmoperpt} #each observation's Hyvarinen divergence weighted by w
  vargradsmo <- cov(gradsmoperpt)

  sensinv <- solve(sens)
  Ginfinv <- sensinv %*% vargradsmo %*% sensinv #inverse of the Godambe information matrix, also called the sandwich information matrix
  out <- Ginfinv/nrow(gradsmoperpt)
  return(out) 
}
