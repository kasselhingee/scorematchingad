#' @title Score matching estimator for the von-Mises Fisher distribution
#' @family directional model estimators
#' @description
#' We recommend using [`movMF::movMF()`] over the function here, `vMF()`, for estimating the von Mises Fisher distribution. 
#' 
#' \insertCite{mardia2016sc;textual}{scorecompdir} suggest score matching for the von Mises Fisher distribution to avoid the difficult-to-compute normalising constant of the von Mises Fisher distributions.
#' \insertCite{mardia2016sc;textual}{scorecompdir} suggested both full score matching and a hybrid estimator, where by the mean direction is estimated via the maximum-likelihood approach and the concentration is estimated via score matching.
#' However, \insertCite{mardia2016sc;textual}{scorecompdir} found, that these score matching estimators were not as a efficient as the maximum likelihood estimators (e.g. [`movMF::movMF()`] ) on the circle.
#' Some of our own experiments suggest this is also the case for higher dimensions, despite the more complex approximations required for the maximum likelihood estimators.

#' @section von Mises Fisher Model: 
#' The von Mises Fisher density is proportional to
#' \deqn{\exp(\kappa \mu^T z),}
#' where \eqn{z} is on a unit sphere,
#' \eqn{\kappa} is termed the *concentration*,
#' and \eqn{\mu} is the *mean direction vector*.
#' The effect of the \eqn{\mu} and \eqn{\kappa} can be decoupled in a sense \insertCite{@p169, @mardia2000di}{scorecompdir}, allowing for estimating \eqn{\mu} and \eqn{\kappa} separately.
#'

#' @details
#' The method "smfull" uses score matching to estimate the vector \eqn{\kappa \mu}.
#' The method "Mardia" uses [`vMF_stdY()`] and [`vMF_kappa()`] to estimate \eqn{\kappa} and \eqn{\mu} seperately.
#' @param Y A matrix of multivariate observations in Cartesian coordinates. Each row is a measurement.
#' @param paramvec `smfull` method only: Optional. A vector of same length as the dimension, representing the elements of the vector that is the element-wise (Hadamard) product \eqn{\kappa \mu}. 
#' @param control Control parameters passed to [`Rcgmin::Rcgmin()`].
#' @param method Either "Mardia" for the hybrid score matching estimator from \insertCite{mardia2016sc;textual}{scorecompdir}.
#'  or "smfull" for the full score matching estimator.
#' @param w An optional vector of weights for each measurement in `Y`
#' @references
#' \insertAllCited{}
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
vMF <- function(Y, paramvec = NULL, method = "smfull", control = default_Rcgmin(), w = rep(1, nrow(Y))){
  fit <- NULL
  if (method == "smfull"){
    if (is.null(paramvec)){
      paramvec <- rep(NA, ncol(Y))
    }
    fit <- vMF_full(Y, paramvec, control = control, w=w)
  }
  if (method %in% c("Mardia")){
    if (!is.null(paramvec)){if (any(!is.na(paramvec))){stop("Mardia estimator cannot fix any elements of paramvec")}}
    fit <- vMF_Mardia(Y, control = control, w=w)
  }
  if (is.null(fit)){stop(sprintf("Method '%s' is not valid", method))}

  return(fit)
}


#for vMF_Mardia startk must be the value of the k parameter
vMF_Mardia <- function(sample, control = default_Rcgmin(), w = rep(1, nrow(sample))){
  mu <- vMF_m(sample, w = w)
  samplestd <- vMF_stdY(sample, m = mu, w = w)
  # check: mustd <- colMeans(samplestd); mustd <- mustd / sqrt(sum(mustd^2))
  kappaest <- vMF_kappa(Y = samplestd, w = w, 
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


vMF_full <- function(sample, usertheta, control = default_Rcgmin(), w = NULL){
  p <- ncol(sample)
  stopifnot(length(usertheta) == p)

  tapes <- buildsmotape("sph","identity", "sph", "vMF",
                        rep(1, p)/sqrt(p), 
                        usertheta = usertheta,
                        weightname = "ones",
                        verbose = FALSE)
  out <- cppad_closed(tapes$smotape, Y = sample, w=w)
  theta <- t_fu2t(out$est, usertheta)

  if (isa(out$SE, "numeric")){
    SE <- t_fu2t(out$SE, 0 * usertheta)
    SE <- list(paramvec = SE)
  } else {SE <- out$SE}
  return(list(
    est = list(paramvec = theta,
               k = sqrt(sum(theta^2)),
               m = theta / sqrt(sum(theta^2))),
    SE = SE,
    info = out
  ))
}

#vMF_paramvec <- function(m, k){
#  return(k*m)
#}

vMF_fromparamvec <- function(paramvec){
  k <- sqrt(sum(paramvec^2))
  m <- paramvec / k
  return(list(m = m, k = k))
}
