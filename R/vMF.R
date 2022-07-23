#' @title Score matching estimate for the von-Mises Fisher distribution
#' @param sample Samples from the vMF distribution in cartesian Rd coordinates. Each row is a measurement.
#' @param km Optional for the `smfull` method. A vector of same length as the dimension, representing the elements of the concentraction * Fisher direction vector.
#' If supplied, the non-NA elements are fixed.
#' @param control Control parameters passed to `Rcgmin::Rcgmin()` and eventually `FixedPoint::FixedPoint()`
#' @param method Either `Mardia` for the hybrid score matching estimate from Mardia et al 2016
#'  or `smfull` for a full score matching estimate.
#'  `Mardia_robustsm` will apply robust score matching to only the `kappa` component of the Mardia estimator.
#' @param w An optional vector of weights for each measurement in `sample`
#' @param cW Optional. If supplied then robust estimation using the Windham weights method is applied, and the value of `cW` is the robustness tuning constant.
#' @examples
#' sample <- Directional::rvmf(100, c(1, -1) / sqrt(2), 3)
#' vMF(sample, method = "smfull")
#' vMF(sample, method = "Mardia")
#' @export
vMF <- function(sample, km = NULL, method = "smfull", control = c(default_Rcgmin(), default_FixedPoint()), w = rep(1, nrow(sample)), cW = NULL){
  firstfit <- NULL
  controls <- splitcontrol(control)
  if (method == "smfull"){
    if (is.null(km)){
      km <- rep(NA, ncol(sample))
    }
    starttheta <- t_u2s_const(km, 0.1)
    isfixed <- t_u2i(km)
    firstfit <- vMF_full(sample, starttheta, isfixed, control = controls$Rcgmin, w=w)
  }
  if (method %in% c("Mardia", "Mardia_robustsm")){
    if (method == "Mardia_robustsm"){stopifnot(!is.null(cW))} #for robust score matching, need to use  cW
    stopifnot(is.null(km))
    firstfit <- vMF_Mardia(sample, startk = 10, control = controls$Rcgmin, w=w)
    isfixed <- rep(FALSE, ncol(sample)) #for Windham robust estimation, if it is used
  }
  if (is.null(firstfit)){stop(sprintf("Method '%s' is not valid", method))}

  if(is.null(cW)){return(firstfit)}

  ###### Extra stuff for robust fit with Windham weights
  inWW <- rep(TRUE, ncol(sample))
  ldenfun <- function(sample, theta){ #here theta is km
    k <- sqrt(sum(theta^2))
    m <- theta/k
    return(drop(Directional::dvmf(sample, k, m, logden = TRUE)))
  }

  if (method == "smfull"){
    estimator <- function(Y, starttheta, isfixed, w){
      km <- rep(NA, ncol(Y))
      km[isfixed] <- starttheta[isfixed]
      out <- vMF_full(Y, starttheta, isfixed,
                      control = controls$Rcgmin, w=w)
      return(out$km)
    }
  } else if (method == "Mardia"){
    estimator <- function(Y, starttheta, isfixed, w){
      startk <- sqrt(sum(starttheta^2))
      out <- vMF_Mardia(Y, startk, control = controls$Rcgmin, w=w)
      return(out$km)
    }
  } else if (method == "Mardia_robustsm"){
    sample <- vMF_stdY(sample, firstfit$m) #standardise sample
    estimator <- function(Y, starttheta, isfixed, w){
      startk <- sqrt(sum(starttheta^2))
      m <- starttheta/startk
      out <- vMF_kappa(Y, startk, control = controls$Rcgmin, w=w)
      return(m * out$k)
    }
  }
  if (is.null(estimator)){stop(sprintf("Method '%s' is not valid", method))}

  est <- windham_raw(prop = sample,
                     cW = cW,
                     ldenfun = ldenfun,
                     estimatorfun = estimator,
                     starttheta = firstfit$km,
                     isfixed = isfixed,
                     inWW = inWW,
                     originalcorrectionmethod = TRUE,
                     fpcontrol = controls$fp)
  est$km <- est$theta
  est$theta <- NULL
  est$k <- sqrt(sum(est$km^2))
  est$m <- est$km/est$k
  return(est)
}

vMF_stdY <- function(Y, m = NULL, w = NULL){
  if(is.null(m)){
    m <- apply(Y, MARGIN = 2, weighted.mean, w)
    m <- mu/sqrt(sum(m^2))
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
  samplestd <- vMF_stdY(sample, m = mu)
  # check: mustd <- colMeans(samplestd); mustd <- mustd / sqrt(sum(mustd^2))
  kappainfo <- vMF_kappa(samplestd, startk, isfixed = isfixed, control = control, w = w)
  return(list(
    k = kappainfo$k,
    m = mu,
    km =  kappainfo$k * mu,
    SE = list(k = kappainfo$SE),
    sminfo = kappainfo
  ))
}

# this function assumes a standardised data set with mean direction equal to the northpole
vMF_kappa <- function(Y, startk, isfixed = FALSE, control = default_Rcgmin(), w = rep(1, nrow(Y))){
  # do estimate, where all but the first component of theta are fixed at zero
  # because kappa * e1 = (kappa, 0, 0, 0, ...)
  if (!isfixed){ #as if k isn't supplied
    p <- ncol(Y)
    tapes <- buildsmotape_internal("Snative", "vMF",
                                   rep(1, p)/sqrt(p), 
                                   starttheta = c(startk, rep(0, p-1)), 
                                   isfixed = c(FALSE, rep(TRUE, p-1)),
                                   weightname = "ones",
                                   verbose = FALSE)
    sminfo <- cppadest(tapes$smotape, startk, Y, control = control, w = w)
    k <- sminfo$par
    SE <- sminfo$SE
  } else {
    sminfo <- NULL
    k <- startk
    SE <- 0
  }
  return(list(
    k = k,
    SE = SE,
    sminfo = sminfo
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
    km = theta,
    k = sqrt(sum(theta^2)),
    m = theta / sqrt(sum(theta^2)),
    SE = list(km = SE),
    sminfo = out
  ))
}
