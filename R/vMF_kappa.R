#' @title Estimate the concentration for a von Mises Fisher distribution
#' @description Using score matching, estimates the concentration \eqn{\kappa} from a sample with mean direction of `c(1, 0, 0, ..., 0)`. 
#' Results by \insertCite{mardia2016sc;textual}{scorecompdir} and some experiments of our own suggest that good implementations of the maximum likelihood estimator (e.g. [`movMF::movMF()`] ) will out perform `vMF_kappa()`.
#'
#' Often a sample with mean direction of `c(1, 0, 0, ...., 0)` is created by estimating the mean direction and rotating the data such that the mean direction equals `c(1, 0, 0, ...)`.
#' Performing this mean direction estimate, then rotation, then estimating concentration with score matching correponds to the hybrid estimator by \insertCite{mardia2016sc;textual}{scorecompdir}.
#' @inherit vMF sections
#' @family Mardia hybrid estimators
#' @references \insertAllCited()
#' @details
#' The function `vMF_kappa()` estimates *only* the concentration \eqn{\kappa}, and assumes that \eqn{\mu} is \eqn{(1, 0, 0, ..., 0)}.
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
    tapes <- buildsmotape("Snative", "vMF",
                          utape = rep(1, p)/sqrt(p),
                          usertheta = c(NA, rep(0, p-1)))
    sminfo <- cppad_closed(tapes$smotape, Y, w = w)
    k <- sminfo$est
    SE <- sminfo$SE
  return(list(
    k = k,
    SE = SE,
    sminfo = sminfo
  ))
}

