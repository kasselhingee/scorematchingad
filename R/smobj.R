#' @title Compute value of score matching objective for multiple observations
#' @description The value of the score matching objective for a given beta
#' @param smofun A tape of the score matching objective calculation
#' @param theta A parameter set
#' @param utabl A matrix of observations, each row being an observation.
#' @export
smobj <- function(smofun, theta, utabl){
  sc_perpt <- lapply(1:nrow(utabl), function(i){
    scobj <- pForward0(smofun, theta, utabl[i,])
    return(scobj)
  })
  scmo <- mean(unlist(sc_perpt))
  return(scmo)
}

#' @describeIn smobj The gradient of the score matching objective function at given beta
#' @export
smobjgrad <- function(smofun, theta, utabl){
  sc_perpt <- lapply(1:nrow(utabl), function(i){
    scobj <- pJacobian(smofun, theta, utabl[i,])
    return(scobj)
  })
  scmo <- colMeans(do.call(rbind, sc_perpt))
  return(scmo)
}

smobjhess <- function(smofun, theta, utabl){
  sc_perpt <- lapply(1:nrow(utabl), function(i){
    scobj <- pHessian(smofun, theta, utabl[i,])
    return(scobj)
  })
  scmo <- colMeans(do.call(rbind, sc_perpt))
  dim(scmo) <- rep(length(theta), 2)
  return(scmo)
}

smestGinf <- function(smofun, est, utabl){ #also called the sandwich information matrix
  warning("Theory for this method not confirmed. Ideas from the Section 2.2 of Varin et al (2011) 'An overview of composite likelihood methods'")
  sens <- -smobjhess(smofun, est, utabl)
  gradsmoperpt <- lapply(1:nrow(utabl), function(i){
    diff <- pJacobian(smofun, est, utabl[i,])
    return(diff)
  })
  vargradsmo <- cov(do.call(rbind, gradsmoperpt)) #SAMPLE estimate of population VARIANCE of gradsmo
  return(sens %*% solve(vargradsmo) %*% sens)
}

smestSE <- function(smofun, est, utabl){
  warning("Theory for this method is not confirmed. Ideas from the Section 2.2 and Section 2.3 of Varin et al (2011) 'An overview of composite likelihood methods'")
  sens <- -smobjhess(smofun, est, utabl)
  gradsmoperpt <- lapply(1:nrow(utabl), function(i){
    diff <- pJacobian(smofun, est, utabl[i,])
    return(diff)
  })
  vargradsmo <- cov(do.call(rbind, gradsmoperpt)) #SAMPLE estimate of population VARIANCE of gradsmo
  sensinv <- solve(sens)
  Ginfinv <- sensinv %*% vargradsmo %*% sensinv #also called the sandwich information matrix
  return(sqrt(Ginfinv / nrow(utabl))) #results now in same units as estimates
}
