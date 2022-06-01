#' @title Find the optimal parameter for the score matching objective
#' @description Using a numerical optimiser and `smobj()` and `smobjgrad()`, find the parameter set with minimal `smobj` value.
#' This function lightly wraps `Rcgmin::Rcgmin()` with some results checking.
#' @param smofun A tape of the score matching objective calculation
#' @param theta The starting parameter set
#' @param utabl A matrix of observations, each row being an observation.
#' @param uboundary Optional. A matrix of observations close to or on the boundary.
#' The contribution of the score matching objective (and corresponding derivatives) of these measurements
#' will be approximated using Taylor approximation from the corresponding `boundaryapprox` matrix.
#' @param boundaryapprox Optional. A matrix of locations in the manifold associated with each measurement in `boundary`.
#' For best results these locations should be further from the manifold boundary and close to their corresponding measurements.
#' Taylor approximation around the rows of `boundaryapprox` will be used to approximate the score matching objective for these measurements.
#' @param control Optional argument passed to `Rcgmin::Rcgmin()`.
#' The default `list(tol = 1E-20)` means the optimisation won't end until the squared size of the gradiant
#'  at the estimate is less than 1E-20.
#' @param approxcentres Optional. A function that takes the measurments that are `TRUE` according to `isboundary`
#' @return The output from `Rcgmin::Rcgmin()`, the squared size of the gradient at the estimate, and the standard error estimates by `smestSE()`.
#' @export
smest <- function(smofun, theta, utabl, control = list(tol = 1E-20), uboundary = NULL, boundaryapprox = NULL){
  if (!(is.null(uboundary) && is.null(boundaryapprox))){
    stopifnot((!is.null(uboundary)) && (!is.null(boundaryapprox))) #both need to be supplied
  }
  out <- Rcgmin::Rcgmin(par = theta,
                        fn = function(theta){smobj(smofun, theta, utabl)},
                        gr = function(theta){smobjgrad(smofun, theta, utabl)},
                        control = control)
  if (out$convergence != 0){warning("Optimisation did not converge.")}
  out$SE <- try({smestSE(smofun, out$par, utabl)})
  out$sqgradsize <- sum(smobjgrad(smofun, out$par, utabl)^2)
  return(out)
}

smest_simplex <- function(smofun, theta, utabl, control = list(tol = 1E-20), shiftsize){
  splittbl <- splitutable_nearbdry(utabl, boundary = "simplex", bdrythreshold = shiftsize, shiftsize = shiftsize)
  smofun_u <- swapDynamic(smofun, rep(1/ncol(utabl), ncol(utabl)), theta * 0 - 0.1) #don't use a boundary point here!
  Jsmofun_u <- pTapeJacobianSwap(smofun, theta * 0 - 0.1, rep(1/ncol(utabl), ncol(utabl)))
  out <- Rcgmin::Rcgmin(par = theta,
                        fn = function(theta){smobj2(smofun, smofun_u, theta,
                                                    splittbl$interior, splittbl$bdry, splittbl$acentres, approxorder = 100)},
                        gr = function(theta){smobjgrad2(smofun, Jsmofun_u, theta,
                                                        splittbl$interior, splittbl$bdry, splittbl$acentres, approxorder = 100)},
                        control = control)
  if (out$convergence != 0){warning("Optimisation did not converge.")}
  return(out)
}



