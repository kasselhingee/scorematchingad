#' @title Score matching estimate for the von-Mises Fisher distribution
#' @details The parameter vector for vMF is the `K * m` where `k` is the concentration, `m` is the mean direction.
#' @param sample Samples from the vMF distribution in cartesian Rd coordinates. Each row is a measurement.
#' @param paramvec Optional for the `smfull` method. A vector of same length as the dimension, representing the elements of the concentraction * Fisher direction vector. It is `k * m` where `k` is the concentration, `m` is the mean direction.
#' If supplied, the non-NA elements are fixed. It is easy to generate `paramvec` using [vMF_paramvec()].
#' @param control Control parameters passed to `Rcgmin::Rcgmin()`.
#' @param method Either `Mardia` for the hybrid score matching estimate from Mardia et al 2016
#'  or `smfull` for a full score matching estimate.
#' @param w An optional vector of weights for each measurement in `sample`
#' @examples
#' set.seed(12342)
#' Y <- movMF::rmovMF(1000, 100 * c(1, 1) / sqrt(2))
#' movMF::movMF(Y, 1) #maximum likelihood estimate
#' vMF(Y, method = "smfull")
#' vMF(Y, method = "Mardia")
#' @return
#' A list of `est`, `SE` and `info`.
#' `est` contains the estimates in vector form, `paramvec`, and with user friendly names `k` and `m`.
#' `SE` contains estimates of the standard errors if computed by the estimating method. Standard error estimates of `k` and `m` are not computed.
#' `info` contains a variety of information about the model fitting procedure.
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
  mu <- apply(sample, MARGIN = 2, weighted.mean, w)
  mu <- mu/sqrt(sum(mu^2))
  samplestd <- vMF_stdY(sample, m = mu, w = w)
  # check: mustd <- colMeans(samplestd); mustd <- mustd / sqrt(sum(mustd^2))
  kappaest <- vMF_kappa(Y = samplestd, w = w, paramvec_start = startk, 
                        control = control)
  return(list(
    est = list(paramvec = kappainfo$k * mu,
               k = kappaest$k,
               m = mu),
    SE = list(paramvec = "Not calculated.",
              k = kappaest$SE,
              m = "Not calculated."),
    info = kappainfo
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
