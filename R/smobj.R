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
