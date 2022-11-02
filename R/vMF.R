#' @title Score matching estimator for the von-Mises Fisher distribution
#' @description
#' We recommend using [`movMF::movMF()`] over the function here, `vMF()`, for estimating the von Mises Fisher distribution.
#' 
#' \insertCite{@mardia2016sc} suggest score matching for the von Mises Fisher distribution to avoid the difficult-to-compute normalising constant of the von Mises Fisher distributions.
#' \insertCite{@mardia2016sc} suggested both full score matching and a hybrid estimator, where by the mean direction is estimated via the maximum-likelihood approach and the concentration is estimated via score matching.
#' However, \insertCite{@mardia2016sc} found, consistent with some of our preliminary experiments, that these score matching estimators were not as a efficient as the maximum likelihood estimators (e.g. [`movMF::movMF()`]).
#' @section von Mises Fisher Model: 
#' The von Mises Fisher density is proportional to
#' \deqn{\exp(\kappa \mu^T z),}
#' where \eqn{z} is on a unit sphere,
#' \eqn{kappa} is termed the *concentration*,
#' and \eqn{\mu} is the *mean direction vector*.
#' The effect of the \eqn{\mu} and \eqn{\kappa} can be decoupled in a sense \insertCite{@mardia2000??}{scorecompdir}, allowing for estimating \eqn{\mu} and \eqn{\kappa} separately.
#'
#' The method "smfull" uses score matching to estimate the vector \eqn{\kappa \mu}.
#' The method "Mardia" uses [`vMF_stdY()`] and [`vMF_kappa()`] to estimate the parameters of the von Mises Fisher distribution.

#' @details
#' The parameter vector for vMF is the `K * m` where `k` is the concentration, `m` is the mean direction.
#' @param Y A matrix of multivariate observations in Cartesian coordinates. Each row is a measurement.
#' @param paramvec `smfull` method only: Optional. A vector of same length as the dimension, representing the elements of the vector \eqn{\kappa \mu}. 
#' If supplied, the non-NA elements are fixed. It is easy to generate `paramvec` using [`vMF_paramvec()`].
#' @param control Control parameters passed to [`Rcgmin::Rcgmin()`].
#' @param method Either "Mardia" for the hybrid score matching estimator from \insertCite{@mardia2016sc}{scorecompdir}.
#'  or "smfull" for the full score matching estimator.
#' @param w An optional vector of weights for each measurement in `Y`
#' @examples
#' set.seed(12342)
#' Y <- movMF::rmovMF(1000, 100 * c(1, 1) / sqrt(2))
#' movMF::movMF(Y, 1) #maximum likelihood estimate
#' vMF(Y, method = "smfull")
#' vMF(Y, method = "Mardia")
#' @return
#' A list of `est`, `SE` and `info`.
#'  * `est` contains the estimates in vector form, `paramvec`, and with user friendly names `k` and `m`.
#'  * `SE` contains estimates of the standard errors if computed.
#'  * `info` contains a variety of information about the model fitting procedure and results.
#' @export
vMF <- function(Y, paramvec = NULL, method = "smfull", control = default_Rcgmin(), w = rep(1, nrow(Y)), paramvec_start = NULL){
  firstfit <- NULL
  if (method == "smfull"){
    if (is.null(paramvec)){
      paramvec <- rep(NA, ncol(Y))
    }
    if (is.null(paramvec_start)){ starttheta <- t_u2s_const(paramvec, 0.1) }
    else {starttheta <- t_us2s(paramvec, paramvec_start) }
    isfixed <- t_u2i(paramvec)
    firstfit <- vMF_full(Y, starttheta = starttheta, isfixed, control = control, w=w)
  }
  if (method %in% c("Mardia")){
    if (!is.null(paramvec)){if (any(!is.na(paramvec))){stop("Mardia estimator cannot fix any elements of paramvec")}}
    if (is.null(paramvec_start)){ startk <- 10 }
    else {startk <- sqrt(sum(paramvec_start^2))}
    firstfit <- vMF_Mardia(Y, startk = startk, control = control, w=w)
    isfixed <- rep(FALSE, ncol(Y)) #for Windham robust estimation, if it is used
  }
  if (is.null(firstfit)){stop(sprintf("Method '%s' is not valid", method))}

  return(firstfit)
}


#for vMF_Mardia startk must be the value of the k parameter
vMF_Mardia <- function(sample, startk, isfixed = FALSE, control = default_Rcgmin(), w = rep(1, nrow(sample))){
  stopifnot(length(startk) == 1)
  stopifnot(length(isfixed) == 1)
  stopifnot(!isfixed)
  mu <- vMF_m(sample, w = w)
  samplestd <- vMF_stdY(sample, m = mu, w = w)
  # check: mustd <- colMeans(samplestd); mustd <- mustd / sqrt(sum(mustd^2))
  kappaest <- vMF_kappa(Y = samplestd, w = w, paramvec_start = startk, 
                        control = control)
  return(list(
    est = list(paramvec = kappaest$k * mu,
               k = kappaest$k,
               m = mu),
    SE = list(paramvec = "Not calculated.",
              k = kappaest$SE,
              m = "Not calculated."),
    info = kappaest
  ))
}


vMF_full <- function(sample, starttheta, isfixed, control = default_Rcgmin(), w = NULL){
  p <- ncol(sample)
  stopifnot(length(starttheta) == p)
  stopifnot(length(starttheta) == length(isfixed))

  tapes <- buildsmotape_internal("Snative", "vMF",
                        rep(1, p)/sqrt(p), starttheta, isfixed,
                        weightname = "ones",
                        verbose = FALSE)
  out <- cppadest(tapes$smotape, t_si2f(starttheta, isfixed), sample, control = control, w=w)
  theta <- t_sfi2u(out$par, starttheta, isfixed)

  SE <- t_sfi2u(out$SE, rep(0, length(starttheta)), isfixed)
  return(list(
    est = list(paramvec = theta,
               k = sqrt(sum(theta^2)),
               m = theta / sqrt(sum(theta^2))),
    SE = list(paramvec = SE),
    info = out
  ))
}

vMF_paramvec <- function(m, k){
  return(k*m)
}

vMF_fromparamvec <- function(paramvec){
  k <- sqrt(sum(paramvec^2))
  m <- paramvec / k
  return(list(m = m, k = k))
}
