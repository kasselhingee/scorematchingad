#' @noRd
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

  # fake xtape and dyntape attributes
  if(!isTRUE(nrow(utabl) > 0)){
    warning("Guessing an interior point for taping")
    p <- ncol(utabl)
    eginteriorpt <- rep(1/p, p)
  } else {
    eginteriorpt <- utabl[1, ]
  }
  attr(smofun, "xtape") <- theta
  attr(smofun, "dyntape") <- eginteriorpt
  out <- cppad_search(smotape = smofun, theta = theta, Y = utabl, Yapproxcentres = NA * utabl,
                      w = w, approxorder = approxorder, control = control)
  return(out)
}

# same parameters as cppad_closed
#' @title Score Matching Estimator for Using Conjugate-Gradient Descent
#' @description 
#' Returns the vector when the gradient of the score matching objective is within tolerance of zero.
#' Also estimates standard errors and covariance.
#' Useful when the score matching objective function is not of quadratic form.
#' For score matching objective functions that are quadratic [`cppad_closed()`] will be usually be more accurate and faster.
#' @param Yapproxcentres A matrix of Taylor approximation centres for rows of Y that require approximation. `NA` for rows that do not require approximation.
#' @param smotape A tape of a score matching objective function.
#' The `smotape`'s independent variables are assumed to be the model parameters to fit
#' and the `smotape`'s dynamic parameter is a (multivariate) measurement.
#' @param Y A matrix of multivariate observations. Each row is an observation.
#' @param w Weights for each observation.
#' @param approxorder The order of Taylor approximation to use.
#' @param control Control parameters passed to [`Rcgmin::Rcgmin()`]
#' @details
#' The score matching objective function and gradient of the score matching function are passed to [`Rcgmin::Rcgmin()`]. 
#' Taylor approximations are performed using [`pTaylorApprox()`] for measurements with non-`NA` rows in `Yapproxcentres`.
#' The call to [`Rcgmin::Rcgmin()`] uses the *sum* of observations (as opposed to the mean) to reduce floating point inaccuracies. This has implications for the meaning of the control parameters passed to `Rcgmin()` (e.g. `tol`). The results are converted into averages so the use of sums can be ignored when not setting control parameters, or studying the behaviour of Rcgmin. 
#'
#' Standard errors are only computed when the weights are constant, and use the Godambe information matrix (aka sandwich method) (*would like a good reference for this here*).
#' The sensitivity matrix \eqn{G} is estimated as
#' the negative of the average over the Hessian of `smotape` evaluated at each observation in `Y`.
# \deqn{\hat{G(\theta)} = \hat{E} -H(smo(\theta;Y))),}
# where \eqn{smo} is the score matching objective function represented by `smotape`,
# \eqn{H} is the Hessian with respect to \eqn{\theta}, which is constant for quadratic-form functions,
# 
#' The variability matrix \eqn{J} is then estimated as
#' the sample covariance (denominator of \eqn{n-1}) of the gradiant of `smotape` evaluated at each of the observations in `Y` for the estimated \eqn{\theta}.
# \deqn{\hat{J}(\theta) = var(grad(w smo(\theta;Y))),}

#' The variance of the estimator is then estimated as
#' \eqn{G^{-1}JG^{-1}/n,}
# \deqn{\hat{G}(\theta)^{-1}\hat{J}(\theta)\hat{G}(\theta)^{-1}/n,}
#' where `n` is the number of observations.
#' @export
cppad_search <- function(smotape, theta, Y, Yapproxcentres = NA * Y, w = rep(1, nrow(Y)), approxorder = 10, control = default_Rcgmin()){
  Jsmofun <- pTapeJacobian(smotape, attr(smotape, "xtape"), attr(smotape, "dyntape"))
  Hsmofun <- pTapeJacobian(Jsmofun, attr(smotape, "xtape"), attr(smotape, "dyntape"))
  
  smofun_u <- swapDynamic(smotape, attr(smotape, "dyntape"), attr(smotape, "xtape")) #don't use a boundary point here!
  Jsmofun_u <- swapDynamic(Jsmofun, attr(smotape, "dyntape"), attr(smotape, "xtape"))

  smoobj <- function(atheta){
    tape_eval_wsum(smofun_u, xmat = Y, pmat = atheta, w = w, xcentres = Yapproxcentres, approxorder = approxorder)
  }
  smograd <- function(atheta){
    tape_eval_wsum(Jsmofun_u, xmat = Y, pmat = atheta, w = w, xcentres = Yapproxcentres, approxorder = approxorder)
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
    sqrt(diag(sme_estvar(smotape, estimate = out$par, Y = Y, Yapproxcentres = Yapproxcentres, approxorder = approxorder)))
  })
  gradatest <- smograd(out$par) / sum(w)
  out$sqgradsize <- sum(gradatest^2)
  return(out)
}


