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

smobj2 <- function(smofun, smofun_u, theta, interior, bdry, acentres, approxorder = 100){
  scmo1 <- 0
  if (nrow(interior) > 0){scmo1 <- smobj(smofun, theta, interior)}

  scmo2 <- 0
  if (nrow(bdry) > 0){
    scmo2 <- smobj_approx(smofun_u, theta, bdry, acentres, approxorder = approxorder)
  }
  scmo <- weighted.mean(c(scmo1, scmo2), w = c(nrow(interior), nrow(bdry)))
  return(scmo)
}

smobjgrad2 <- function(smofun, Jsmofun_u, theta, interior, bdry, acentres, approxorder = 100){
  scmo1 <- 0
  if (nrow(interior) > 0){scmo1 <- smobjgrad(smofun, theta, interior)}

  scmo2 <- 0
  if (nrow(bdry) > 0){
    scmo2 <- smobjgrad_approx(Jsmofun_u, theta, bdry, acentres, approxorder = approxorder)
  }
  scmo <- weighted.mean(c(scmo1, scmo2), w = c(nrow(interior), nrow(bdry)))
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


smobj_approx <- function(smofun_u, theta, utabl, acentres, approxorder = 100){
  stopifnot(nrow(utabl) == nrow(acentres))
  sc_perpt <- lapply(1:nrow(utabl), function(i){
    scobj <- pTaylorApprox(smofun_u, utabl[i,], acentres[i, ], theta, order = approxorder)
    return(scobj)
  })
  scmo <- mean(unlist(sc_perpt))
  return(scmo)
}

smobjgrad_approx <- function(Jsmofun_u, theta, utabl, acentres, approxorder = 100){
  stopifnot(nrow(utabl) == nrow(acentres))
  sc_perpt <- lapply(1:nrow(utabl), function(i){
    scobj <- pTaylorApprox(Jsmofun_u, utabl[i,], acentres[i, ], theta, order = approxorder)
    return(scobj)
  })
  scmo <- colMeans(unlist(sc_perpt))
  return(scmo)
}
