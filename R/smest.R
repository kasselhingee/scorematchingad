#' @noRd
#' @title Find the optimal parameter for the score matching objective
#' @description Using a numerical optimiser and `smobj()` and `smobjgrad()`, find the parameter set with minimal `smobj` value.
#' This function lightly wraps `Rcgmin::Rcgmin()` with some results checking.
#' @details The call to `Rcgmin()` uses the *sum* of observations (as opposed to the mean) to reduce floating point inaccuracies. This has implications for the meaning of the control parameters passed to `Rcgmin()` (e.g. `tol`). The results are converted into averages as appropriate so the use of sums (as opposed to averages) can be ignored when not setting control parameters, or studying the behaviour of Rcgmin. 
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
cppadest <- function(smofun, theta, utabl, control = default_Rcgmin(), uboundary = NULL, boundaryapprox = NULL, approxorder = NULL, w = NULL, wboundary = NULL){
  if (!(is.null(uboundary) && is.null(boundaryapprox))){
    stopifnot((!is.null(uboundary)) && (!is.null(boundaryapprox))) #both need to be supplied
    stopifnot(nrow(uboundary) == nrow(boundaryapprox))
    stopifnot(ncol(uboundary) == ncol(boundaryapprox))
    stopifnot(ncol(uboundary) == ncol(utabl))
  }
  
  # convert old-style boundary specification into new style
  boundaryapprox <- rbind(NA * utabl, boundaryapprox)
  utabl <- rbind(utabl, uboundary)
  w <- c(w, wboundary)
  if (is.null(w)){w <- rep(1, nrow(utabl))}

  # generate tapes with respect to the measurement (parameter is dynamic)
  if(!isTRUE(nrow(utabl) > 0)){
    warning("Guessing an interior point for taping")
    p <- ncol(utabl)
    eginteriorpt <- rep(1/p, p)
  } else {
    eginteriorpt <- utabl[1, ]
  }
  Jsmofun <- pTapeJacobian(smofun, theta, eginteriorpt)
  Hsmofun <- pTapeJacobian(Jsmofun, theta, eginteriorpt)
  
  smofun_u <- swapDynamic(smofun, eginteriorpt, theta) #don't use a boundary point here!
  Jsmofun_u <- swapDynamic(Jsmofun, eginteriorpt, theta)
  Hsmofun_u <- swapDynamic(Hsmofun, eginteriorpt, theta)

  smoobj <- function(atheta){
    tape_eval_wsum(smofun_u, xmat = utabl, pmat = atheta, w = w, xcentres = boundaryapprox, approxorder = approxorder)
  }
  smograd <- function(atheta){
    tape_eval_wsum(Jsmofun_u, xmat = utabl, pmat = atheta, w = w, xcentres = boundaryapprox, approxorder = approxorder)
  }

  out <- Rcgmin::Rcgmin(par = theta,
                        fn = smoobj,
                        gr = smograd,
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
  out$value <- out$value / sum(w)

  out$SE <- try({
    sqrt(diag(sme_estvar(smofun, estimate = out$par, Y = utabl, Yapproxcentres = boundaryapprox, approxorder = approxorder)))
  })
  gradatest <- smograd(out$par)
  out$sqgradsize <- sum(gradatest^2)
  return(out)
}




