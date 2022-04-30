#' @title Find the optimal parameter for the score matching objective
#' @description Using a numerical optimiser and `smobj()` and `smobjgrad()`, find the parameter set with minimal `smobj` value.
#' This function lightly wraps `Rcgmin::Rcgmin()` with some results checking.
#' @param smofun A tape of the score matching objective calculation
#' @param theta The starting parameter set
#' @param utabl A matrix of observations, each row being an observation.
#' @param control Optional argument passed to `Rcgmin::Rcgmin()`.
#' The default `list(tol = 1E-20)` means the optimisation won't end until the squared size of the gradiant
#'  at the estimate is less than 1E-20.
#' @return The output from `Rcgmin::Rcgmin()`, the squared size of the gradient at the estimate, and the standard error estimates by `smestSE()`.
#' @export
smest <- function(smofun, theta, utabl, control = list(tol = 1E-20)){
  out <- Rcgmin::Rcgmin(par = theta,
                        fn = function(theta){smobj(smofun, theta, utabl)},
                        gr = function(theta){smobjgrad(smofun, theta, utabl)},
                        control = control)
  if (out$convergence != 0){warning("Optimisation did not converge.")}
  out$SE <- smestSE(smofun, out$par, utabl)
  out$gradsize <- sum(smobjgrad(smofun, out$par, utabl)^2)
  return(out)
}
