#' @title Estimate the concentration for a von Mises Fisher distribution
#' @description Estimates the concentration `k` from a sample with mean direction of `c(1, 0, 0, ...)`. Often such a sample is created by estimating the mean direction and rotating the data such that the mean direction equals `c(1, 0, 0, ...)`.
#' Performing this mean direction estimate, rotation, then estimating concentration with score matching correponds to Mardia et al's \insertCite{@mardia2016sc}{scorecompdir} hybrid estimator for von Mises-Fisher distribution.
#' @details
#' The von Mises Fisher distribution has two parameters, the concentration \eqn{\kappa} and mean direction \eqn{\mu} with density \eqn{f(z; \kappa, \mu)} proportional to
#' \deqn{\exp(\kappa \mu^T z).}
#'  
#' The concentration and mean direction have well defined and in-a-sense independent impact on the distribution \insertCite{p169 @mardia2000di}, which allows for them to be estimated seperately. The function here estimates *only* the concentration, assuming that \eqn{\mu} is \eqn{(1, 0, 0, ..., 0)}.
#' @return
#' A list of `est`, `SE` and `info`.
#' `est` contains the estimate of the concentration in the slot `k` for easy use and `paramvec` for compatibility with other functions in this package.
#' `SE` contains estimates of the standard errors if computed by the estimating method.
#' `info` contains a variety of information about the model fitting procedure.
#' @param paramvec_start The starting value of concentration parameter `k` or kappa for iterative minimisation of the score-matching objective.
#' @param Y A data matrix, each row is an observation.
#' @param w Weights corresponding to each row of `Y`.
#' @param control A list of control arguments passed to [Rcgmin::Rcgmin()].
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

