#' @title Draft standard error estimates for a CppAD-based estimator
#' @description The standard errors of an estimator via the Godambe information matrix (or sandwich information matrix).
#' WARNING: Theory for this method in this situation has not been confirmed by KLH. Ideas from the Section 2.2 and Section 2.3 of Varin et al (2011) 'An overview of composite likelihood methods.
#' @param smofun A tape of the score matching objective calculation
#' @param est The parameter set
#' @param utabl A matrix of observations, each row being an observation.
#' @details It seems likely that `sqrt(n) * (est - true)` will have an asympotically normal distribution with zero mean and a covariance of
#' given by the inverse of the Godambe information matrix, `invG`.
#' @return The square root of the diagonal of the matrix `invG/n`, which is an estimate of estimator SE.
#' @export
smestSE <- function(smofun, theta, utabl, ...){
  sqrt(diag(smestSEsq(smofun, theta, utabl, ...)))
}

#' @describeIn smestSE Returns the matrix `invG/n`.
smestSEsq <- function(smofun, theta, utabl,
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

  if (!is.null(w)){
     if (!all(w == rep(1, length(w)))){
       stop("SE for weighted data not available")
    }
  }

  vargradsmo <- cov(do.call(rbind, grasmoperpt))

  sensinv <- solve(sens)
  Ginfinv <- sensinv %*% vargradsmo %*% sensinv #inverse of the Godambe information matrix, also called the sandwich information matrix
  if (is.null(w)){out <- Ginfinv/length(gradsmoperpt)}
  else {out <- Ginfinv/sum(w)}
  return(out) #results now in same units as estimates
}
