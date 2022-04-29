#' @title Draft standard error estimates for a CppAD-based estimator
#' @description The standard errors of an estimator via the Godambe information matrix (or sandwich information matrix).
#' @param smofun A tape of the score matching objective calculation
#' @param est The parameter set
#' @param utabl A matrix of observations, each row being an observation.
#' @details It seems likely that `sqrt(n) * (est - true)` will have an asympotically normal distribution with zero mean and a covariance of
#' given by the inverse of the Godambe information matrix, `invG`.
#' @return The square root of the diagonal of the matrix `invG/n`, which is an estimate of estimator SE.
#' @export
smestSE <- function(smofun, est, utabl){
  sqrt(diag(smestSEsq(smofun, est, utabl)))
}

#' @describeIn smestSE Returns the matrix `invG/n`.
smestSEsq <- function(smofun, est, utabl){
  warning("Theory for this method is not confirmed. Ideas from the Section 2.2 and Section 2.3 of Varin et al (2011) 'An overview of composite likelihood methods'")
  sens <- -smobjhess(smofun, est, utabl)
  gradsmoperpt <- lapply(1:nrow(utabl), function(i){
    diff <- pJacobian(smofun, est, utabl[i,])
    return(diff)
  })
  vargradsmo <- cov(do.call(rbind, gradsmoperpt)) #SAMPLE estimate of population VARIANCE of gradsmo
  sensinv <- solve(sens)
  Ginfinv <- sensinv %*% vargradsmo %*% sensinv #inverse of the Godambe information matrix, also called the sandwich information matrix
  return(Ginfinv / nrow(utabl)) #results now in same units as estimates
}

