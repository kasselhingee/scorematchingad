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
  out$sqgradsize <- sum(smobjgrad(smofun, out$par, utabl)^2)
  return(out)
}

smest_simplex <- function(smofun, theta, utabl, control = list(tol = 1E-20), shiftsize){
  splittbl <- splitutable_nearbdry(utabl, boundary = "simplex", bdrythreshold = shiftsize, shiftsize = shiftsize)
  smofun_u <- swapDynamic(smofun, rep(1/ncol(utabl), ncol(utabl)), theta * 0 - 0.1) #don't use a boundary point here!
  out <- Rcgmin::Rcgmin(par = theta,
                        fn = function(theta){smobj2(smofun, smofun_u, theta,
                                                    splittbl$interior, splittbl$bdry, splittbl$acentres, approxorder = 100)},
                        # gr = function(theta){smobjgrad(smofun, theta, utabl)},
                        control = control)
  if (out$convergence != 0){warning("Optimisation did not converge.")}
  return(out)
}


splitutable_nearbdry <- function(utabl, boundary = "simplex", bdrythreshold = 1E-15, shiftsize = bdrythreshold){
  onbdry <- rep(FALSE, nrow(utabl))
  if (!is.null(boundary)){
    if (boundary == "simplex"){
      onbdry <- apply(utabl, MARGIN = 1, min) < bdrythreshold
      acentres <- approxcentre(utabl[onbdry, ], shiftsize = shiftsize)
    }
  }
  return(list(
    interior = utabl[!onbdry, ],
    bdry = utabl[onbdry, ],
    acentres = acentres))
}

