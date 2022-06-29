#' @title Score matching estimate for the von-Mises Fisher distribution
#' @param sample Samples from the vMF distribution in cartesian Rd coordinates. Each row is a measurement.
#' @param km Optional for the `smfull` method. A vector of same length as the dimension, representing the elements of the concentraction * Fisher direction vector.
#' If supplied, the non-NA elements are fixed.
#' @param control Control parameters passed to `Rcgmin::Rcgmin()` and eventually `FixedPoint::FixedPoint()`
#' @param method Either `Mardia` for the hybrid score matching estimate from Mardia et al 2016 or `smfull` for a full score matching estimate.
#' @param w An optional vector of weights for each measurement in `sample`
#' @param cW Optional. If supplied then robust estimation using the Windham weights method is applied, and the value of `cW` is the robustness tuning constant.
#' @examples
#' sample <- Directional::rvmf(100, c(1, -1) / sqrt(2), 3)
#' vMF(sample, method = "smfull")
#' vMF(sample, method = "Mardia")
#' @export
vMF <- function(sample, km = NULL, method = "smfull", control = list(tol = 1E-20), w = rep(1, nrow(sample)), cW = NULL){
  out <- NULL
  if (method == "smfull"){
    if (is.null(km)){
      km <- rep(NA, ncol(sample))
    }
    starttheta <- t_u2s_const(km, 0.1)
    isfixed <- t_u2i(km)
    firstfit <- vMF_full(sample, starttheta, isfixed, control = control, w=w)
  }
  if (method == "Mardia"){
    stopifnot(is.null(km))
    firstfit <- vMF_Mardia(sample, startk = 0.1, control = control, w=w)
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
      out <- vMF_full(sample, starttheta, isfixed, 
                      control = control, w=w)
      return(out$km)
    }
  } else if (method == "Mardia"){
    estimator <- function(Y, starttheta, isfixed, w){
      startk <- sqrt(sum(starttheta^2))
      out <- vMF_Mardia(sample, startk, control = control, w=w)
      return(out$km)
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
                     originalcorrectionmethod = TRUE)
  return(est)
}


#for vMF_Mardia startk must be the value of the k parameter
vMF_Mardia <- function(sample, startk, isfixed = FALSE, control = list(tol = 1E-20), w = rep(1, nrow(sample))){
  stopifnot(length(startk) == 1)
  stopifnot(length(isfixed) == 1)

  mu <- apply(sample, MARGIN = 2, weighted.mean, w)
  mu <- mu/sqrt(sum(mu^2))
  Rtrans <- Directional::rotation(mu, c(1, rep(0, length(mu) - 1)))
  samplestd <- sample %*% t(Rtrans)
  # check: mustd <- colMeans(samplestd); mustd <- mustd / sqrt(sum(mustd^2))

  # do estimate, where all but the first component of theta are fixed at zero
  # because kappa * e1 = (kappa, 0, 0, 0, ...)
  if (!isfixed){ #as if k isn't supplied
    p <- ncol(sample)
    tapes <- buildsmotape_internal("Snative", "vMF",
                        rep(1, p)/sqrt(p), startk, isfixed,
                        weightname = "ones",
                        verbose = FALSE)
    sminfo <- smest(tapes$smotape, t_si2f(startk, isfixed), samplestd, control = control, w = w)
    k <- sminfo$par
    SE <- list(k = sminfo$SE)
  } else {
    sminfo <- NULL
    k <- startk
    SE <- list(k = 0)
  }
  return(list(
    k = k,
    m = mu,
    km = k * mu,
    SE = SE,
    sminfo = sminfo
  ))
}

vMF_full <- function(sample, starttheta, isfixed, control = list(tol = 1E-20), w = NULL){
  p <- ncol(sample)
  stopifnot(length(starttheta) == p)
  stopifnot(length(starttheta) == length(isfixed))

  tapes <- buildsmotape_internal("Snative", "vMF",
                        rep(1, p)/sqrt(p), starttheta, isfixed,
                        weightname = "ones",
                        verbose = FALSE)
  out <- smest(tapes$smotape, t_si2f(starttheta, isfixed), sample, control = control, w=w)
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
