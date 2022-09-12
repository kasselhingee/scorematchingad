# @title Find the optimal parameter for the score matching objective
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
#' @return The output from `Rcgmin::Rcgmin()`, the squared size of the gradient at the estimate, and the standard error estimates by `cppadSE()`.
# @export
cppadest <- function(smofun, theta, utabl, control = default_Rcgmin(), uboundary = NULL, boundaryapprox = NULL, approxorder = NULL, w = NULL){
  if (!(is.null(uboundary) && is.null(boundaryapprox))){
    stopifnot((!is.null(uboundary)) && (!is.null(boundaryapprox))) #both need to be supplied
    stopifnot(nrow(uboundary) == nrow(boundaryapprox))
    stopifnot(ncol(uboundary) == ncol(boundaryapprox))
    stopifnot(ncol(uboundary) == ncol(utabl))
    if(!isTRUE(nrow(utabl) > 0)){
      warning("Guessing an interior point for taping")
      p <- ncol(utabl)
      eginteriorpt <- rep(1/p, p)
    } else {
      eginteriorpt <- utabl[1, ]
    }
    smofun_u <- swapDynamic(smofun, eginteriorpt, theta) #don't use a boundary point here!
    Jsmofun_u <- pTapeJacobianSwap(smofun, theta, eginteriorpt)
    Hsmofun_u <- pTapeHessianSwap(smofun, theta, eginteriorpt)
  } else {
    smofun_u <- NULL
    Jsmofun_u <- NULL
    Hsmofun_u <- NULL
  }

  out <- Rcgmin::Rcgmin(par = theta,
                        fn = smobj_sum_b,
                        gr = smobjgrad_sum_b,
                        # function(theta, ...){smobj(smofun, theta, utabl, smofun_u, ...)},
                        # gr = function(theta, ...){smobjgrad(smofun, theta, utabl, Jsmofun_u, ...)},
                        smofun = smofun,
                        utabl = utabl,
                        smofun_u = smofun_u,
                        Jsmofun_u = Jsmofun_u,
                        uboundary = uboundary, boundaryapprox = boundaryapprox,
                        approxorder = approxorder,
                        w = w,
                        stopifnan = TRUE,
                        control = control)
  if (out$convergence == 2){
    if (grepl("Initial point", out$message)){
      stop(paste(out$message, "Initial point was", paste(theta, collapse = " ")))
    } else {
      stop(paste(out$message, "Perhaps smobj() generates a NaN?"))
    }
  }
  if (out$convergence != 0){warning("Optimisation did not converge.")}

  # return results as if averages, not sums were used
  attr(out$par, "normaliser") <- NULL
  out$value <- out$value / attr(out$value, "normaliser")
  attr(out$value, "normaliser") <- NULL


  out$SE <- try({
    cppadSE(
      smofun, theta = out$par, utabl,
      Jsmofun_u = Jsmofun_u, Hsmofun_u = Hsmofun_u,
      uboundary = uboundary, boundaryapprox = boundaryapprox,
      approxorder = approxorder,
      w = w)})
  gradatest <- smobjgrad(smofun, out$par, utabl,
                         Jsmofun_u = Jsmofun_u,
                         uboundary = uboundary, boundaryapprox = boundaryapprox,
                         approxorder = approxorder,
                         w = w)
  out$sqgradsize <- sum(gradatest)^2
  return(out)
}




