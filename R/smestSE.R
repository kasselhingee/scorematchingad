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
smestSE <- function(smofun, est, utabl){
  sqrt(diag(smestSEsq(smofun, est, utabl)))
}

#' @describeIn smestSE Returns the matrix `invG/n`.
smestSEsq <- function(smofun, est, utabl){
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

smestSE_simplex <- function(smofun, est, utabl, shiftsize, approxorder = 100){
  sqrt(diag(smestSEsq_simplex(smofun, est, utabl, shiftsize, approxorder = approxorder)))
}

smestSEsq_simplex <- function(smofun, est, utabl, shiftsize, approxorder = 100){
  splittbl <- splitutable_nearbdry(utabl, boundary = "simplex", bdrythreshold = shiftsize, shiftsize = shiftsize)
  Jsmofun_u <- pTapeJacobianSwap(smofun, est * 0 - 0.1, rep(1/ncol(utabl), ncol(utabl)))
  Hsmofun_u <- pTapeHessianSwap(smofun, est * 0 - 0.1, rep(1/ncol(utabl), ncol(utabl)))
  sens <- -smobjhess2(smofun, Hsmofun_u, est,  splittbl$interior, splittbl$bdry, splittbl$acentres, approxorder = approxorder)

  if (nrow(splittbl$interior) > 0){
    gradsmoperpt_interior <- lapply(1:nrow(splittbl$interior), function(i){
      scobj <- pJacobian(smofun, est, splittbl$interior[i,])
      return(scobj)
    })
  } else {gradsmoperpt_interior <- NULL}

  if (nrow(splittbl$bdry) > 0){
    gradsmoperpt_bdry <- lapply(1:nrow(splittbl$bdry), function(i){
      scobj <- pTaylorApprox(Jsmofun_u, splittbl$bdry[i,], splittbl$acentres[i, ], est, order = approxorder)
      return(scobj)
    })
  } else {gradsmoperpt_bdry <- NULL}
  gradsmoperpt <- c(gradsmoperpt_interior, gradsmoperpt_bdry)
  vargradsmo <- cov(do.call(rbind, gradsmoperpt)) #SAMPLE estimate of population VARIANCE of gradsmo
  sensinv <- solve(sens)
  Ginfinv <- sensinv %*% vargradsmo %*% sensinv #inverse of the Godambe information matrix, also called the sandwich information matrix
  return(Ginfinv / nrow(utabl)) #results now in same units as estimates
}

