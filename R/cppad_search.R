# same parameters as cppad_closed
#' @title Iterative Score Matching Estimator Using Conjugate-Gradient Descent
#' @family generic score matching tools
#' @description 
#' Uses conjugate gradient descent to search for a vector of parameters such that gradient of the score matching objective is within tolerance of zero.
#' Also estimates standard errors and covariance.
#' Useful when the score matching objective function is not of quadratic form;
#' for score matching objective functions that are quadratic [`cppad_closed()`] will be usually be more accurate and faster.
#' @inheritParams cppad_closed
#' @param theta The starting parameter set
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
  Jsmofun <- tapeJacobian(smotape)
  Hsmofun <- tapeJacobian(Jsmofun)
  
  smofun_u <- tapeSwap(smotape) #don't use a boundary point for taping here!
  Jsmofun_u <- tapeSwap(Jsmofun)

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
      stop(paste(out$message, "Perhaps tape_eval_wsum() generates a non-number?"))
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

