#' @title Estimate the concentration for a von Mises Fisher distribution
#' @description Estimates the concentration `k` from a sample with mean direction of `c(1, 0, 0, ...)`. Often such a sample is created by estimating the mean direction and rotating the data such that the mean direction equals `c(1, 0, 0, ...)`.
#' Performing this mean direction estimate, rotation, then estimating concentration with score matching correponds to Mardia et al's hybrid estimators [ref].
#' @details
#' Theoretically, the concentration `k` and mean direction `m` of a von Mises Fisher distribution can be estimated separately [ref Mardia and Judd].
#' ... insert vMF density ...
#' @param paramvec
#' @return
#' A list of `est`, `SE` and `info`.
#' `est` contains the estimate of `k` in `k` for easy use and `paramvec` for compatibility with other functions in this package.
#' `SE` contains estimates of the standard errors if computed by the estimating method.
#' `info` contains a variety of information about the model fitting procedure.
#' @param paramvec_start The starting value of concentration parameter `k` or kappa for iterative minimisation of the score-matching objective.
#' @export
vMF_kappa <- function(Y, w = rep(1, nrow(Y)), paramvec_start = 10, control = default_Rcgmin()){
  stopifnot(length(paramvec_start) == 1)
  # do estimate, where all but the first component of theta are fixed at zero
  # because kappa * e1 = (kappa, 0, 0, 0, ...)
    p <- ncol(Y)
    tapes <- buildsmotape_internal("Snative", "vMF",
                                   rep(1, p)/sqrt(p),
                                   starttheta = c(paramvec_start, rep(0, p-1)),
                                   isfixed = c(FALSE, rep(TRUE, p-1)),
                                   weightname = "ones",
                                   verbose = FALSE)
    sminfo <- cppadest(tapes$smotape, paramvec_start, Y, control = control, w = w)
    k <- sminfo$par
    SE <- sminfo$SE
  return(list(
    k = k,
    SE = SE,
    sminfo = sminfo
  ))
}

