#' @title Score matching estimate for the von-Mises Fisher distribution
#' @details The parameter vector for vMF is the `K * m` where `k` is the concentration, `m` is the mean direction.
#' @param sample Samples from the vMF distribution in cartesian Rd coordinates. Each row is a measurement.
#' @param paramvec Optional for the `smfull` method. A vector of same length as the dimension, representing the elements of the concentraction * Fisher direction vector. It is `k * m` where `k` is the concentration, `m` is the mean direction.
#' If supplied, the non-NA elements are fixed. It is easy to generate `paramvec` using [vMF_paramvec()].
#' @param control Control parameters passed to `Rcgmin::Rcgmin()` and eventually `FixedPoint::FixedPoint()`
#' @param method Either `Mardia` for the hybrid score matching estimate from Mardia et al 2016
#'  or `smfull` for a full score matching estimate.
#'  `Mardia_robustsm` will apply robust score matching to only the `kappa` component of the Mardia estimator.
#' @param w An optional vector of weights for each measurement in `sample`
#' @param cW Optional. If supplied then robust estimation using the Windham weights method is applied, and the value of `cW` is the robustness tuning constant.
#' @examples
#' set.seed(12342)
#' Y <- movMF::rmovMF(1000, 100 * c(1, 1) / sqrt(2))
#' movMF::movMF(Y, 1) #maximum likelihood estimate
#' vMF(Y, method = "smfull")
#' vMF(Y, method = "Mardia")
#' est <- vMF(Y, method = "Mardia_robustsm", cW = 1E-5, control = list(tol = 1E-10, maxit = 100, MaxIter = 100))
#' @return
#' A list of `est`, `SE` and `info`.
#' `est` contains the estimates in vector form, `paramvec`, and with user friendly names `k` and `m`.
#' `SE` contains estimates of the standard errors if computed by the estimating method. Standard error estimates of `k` and `m` are not computed.
#' `info` contains a variety of information about the model fitting procedure.
#' @export
vMF <- function(Y, paramvec = NULL, method = "smfull", control = c(default_Rcgmin(), default_FixedPoint()), w = rep(1, nrow(Y)), cW = NULL){
  firstfit <- NULL
  controls <- splitcontrol(control)
  if (method == "smfull"){
    if (is.null(paramvec)){
      paramvec <- rep(NA, ncol(Y))
    }
    starttheta <- t_u2s_const(paramvec, 0.1)
    isfixed <- t_u2i(paramvec)
    firstfit <- vMF_full(Y, starttheta, isfixed, control = controls$Rcgmin, w=w)
  }
  if (method %in% c("Mardia", "Mardia_robustsm")){
    if (method == "Mardia_robustsm"){stopifnot(!is.null(cW))} #for robust score matching, need to use  cW
    stopifnot(is.null(paramvec))
    firstfit <- vMF_Mardia(Y, startk = 10, control = controls$Rcgmin, w=w)
    isfixed <- rep(FALSE, ncol(Y)) #for Windham robust estimation, if it is used
  }
  if (is.null(firstfit)){stop(sprintf("Method '%s' is not valid", method))}

  if(is.null(cW)){return(firstfit)}

  ###### Extra stuff for robust fit with Windham weights
  ldenfun <- function(Y, theta){ #here theta is km
    k <- sqrt(sum(theta^2))
    m <- theta/k
    return(drop(Directional::dvmf(Y, k, m, logden = TRUE)))
  }

  if (method == "smfull"){
    estimator <- function(Y, starttheta, isfixed, w){
      paramvec <- rep(NA, ncol(Y))
      paramvec[isfixed] <- starttheta[isfixed]
      out <- vMF_full(Y, starttheta, isfixed,
                      control = controls$Rcgmin, w=w)
      return(out$est$paramvec)
    }
  est <- windham_raw(prop = Y,
                     cW = cW * !isfixed,
                     ldenfun = ldenfun,
                     estimatorfun = estimator,
                     starttheta = firstfit$est$paramvec,
                     isfixed = isfixed,
                     originalcorrectionmethod = TRUE,
                     fpcontrol = controls$fp)
  est$paramvec <- est$theta
  est$theta <- NULL
  est$k <- sqrt(sum(est$paramvec^2))
  est$m <- est$paramvec/est$k
  } else if (method == "Mardia"){
    estimator <- function(Y, starttheta, isfixed, w){
      startk <- sqrt(sum(starttheta^2))
      out <- vMF_Mardia(Y, startk, control = controls$Rcgmin, w=w)
      return(out$est$paramvec)
    }
  est <- windham_raw(prop = Y,
                     cW = cW * !isfixed,
                     ldenfun = ldenfun,
                     estimatorfun = estimator,
                     starttheta = firstfit$est$paramvec,
                     isfixed = isfixed,
                     originalcorrectionmethod = TRUE,
                     fpcontrol = controls$fp)
  est$paramvec <- est$theta
  est$theta <- NULL
  est$k <- sqrt(sum(est$paramvec^2))
  est$m <- est$paramvec/est$k
  } else if (method == "Mardia_robustsm"){
    Y <- vMF_stdY(Y, firstfit$est$m, w = w) #standardise sample
    # after standardisation there is only one parameter to fit in the vMF model
    ldenfun <- function(Y, theta){ #here theta is just k
      k <- theta
      m <- c(1, rep(0, ncol(Y) - 1))
      return(drop(Directional::dvmf(Y, k, m, logden = TRUE)))
    }
    estimator <- function(Y, starttheta, isfixed, w){ #starttheta is k in this case
      out <- vMF_kappa_coarse(Y, starttheta, control = controls$Rcgmin, w=w)
      return(out$k)
    }
  est <- windham_raw(prop = Y,
                     cW = cW, #plugs into ldenfun
                     ldenfun = ldenfun,
                     estimatorfun = estimator,
                     starttheta = firstfit$est$k,
                     isfixed = c(FALSE),
                     originalcorrectionmethod = TRUE,
                     fpcontrol = controls$fp)
  est$k <- est$theta
  est$theta <- NULL
  est$m <- firstfit$est$m
  }
  if (is.null(estimator)){stop(sprintf("Method '%s' is not valid", method))}

  return(est)
}

vMF_stdY <- function(Y, m = NULL, w = NULL){
  if(is.null(m)){
    m <- apply(Y, MARGIN = 2, weighted.mean, w)
    m <- m/sqrt(sum(m^2))
  }
  Rtrans <- Directional::rotation(m, c(1, rep(0, length(m) - 1)))
  out <- Y %*% t(Rtrans)
  return(out)
}

#for vMF_Mardia startk must be the value of the k parameter
vMF_Mardia <- function(sample, startk, isfixed = FALSE, control = default_Rcgmin(), w = rep(1, nrow(sample))){
  stopifnot(length(startk) == 1)
  stopifnot(length(isfixed) == 1)
  mu <- apply(sample, MARGIN = 2, weighted.mean, w)
  mu <- mu/sqrt(sum(mu^2))
  samplestd <- vMF_stdY(sample, m = mu, w = w)
  # check: mustd <- colMeans(samplestd); mustd <- mustd / sqrt(sum(mustd^2))
  kappainfo <- vMF_kappa_coarse(samplestd, startk, isfixed = isfixed, control = control, w = w)
  return(list(
    est = list(paramvec = kappainfo$k * mu,
               k = kappainfo$k,
               m = mu),
    SE = list(paramvec = "Not calculated.",
              k = kappainfo$SE,
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
  m <- paramvec / m
  return(list(m = m, k = k))
}
