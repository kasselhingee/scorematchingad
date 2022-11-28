# same parameters as cppad_closed
#' @title Iterative Score Matching Estimator Using Conjugate-Gradient Descent
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
#' @param theta The starting parameter set
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

  # useful to debugging as Rcgmin hides the error codes
  #  and relatively cheap:
  # evaluating above functions at the start point
  stopifnot(is.finite(smoobj(theta)))
  stopifnot(all(is.finite(smograd(theta))))
  # now do the search 
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

  out$SE <- "Not calculated."
  if (isTRUE(all(w[[1]] == w))){
    out$SE <- sqrt(diag(sme_estvar(smotape, estimate = out$par, Y = Y, Yapproxcentres = Yapproxcentres, approxorder = approxorder)))
  }
  gradatest <- smograd(out$par) / sum(w)
  out$sqgradsize <- sum(gradatest^2)
  out$est <- out$par
  out$par <- NULL
  return(out)
}
#' @return A list of
#' + `est` The estimate
#' + `sqgradsize` The squared size of the gradient at the estimate
#' + `SE` The standard error of the estimates when `w` is constant. If `w` changes between observations then this slot contains the string "Not calculated."
#' + additional output from [`Rcgmin::Rcgmin()`]

