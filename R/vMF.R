#' @title Score matching estimate for the von-Mises Fisher distribution
#' @param sample Samples from the vMF distribution in cartesian Rd coordinates. Each row is a measurement.
#' @param km Optional for the `smfull` method. A vector of same length as the dimension, representing the elements of the concentraction * Fisher direction vector.
#' If supplied, the non-NA elements are fixed.
#' @param control Control parameters passed to `Rcgmin::Rcgmin()`
#' @param method Either `Mardia` for the hybrid score matching estimate from Mardia et al 2016 or `smfull` for a full score matching estimate.
#' @examples
#' sample <- movMF::rmovMF(100, 3 * c(1, -1) / sqrt(2))
#' vMF(sample, "smfull")
#' vMF(sample, "Mardia")
#' @export
vMF <- function(sample, km = NULL, method = "smfull", control = list(tol = 1E-20), w = rep(1, nrow(sample))){
  out <- NULL
  if (method == "smfull"){out <- vMF_full(sample, km = km, control = control, w=w)}
  if (method == "Mardia"){
    stopifnot(is.null(km))
    out <- vMF_Mardia(sample, control = control, w=w)
  }
  if (is.null(out)){stop(sprintf("Method '%s' is not valid", method))}
  return(out)
}

vMF_robust <- function(sample, km = NULL, method = "smfull", control = list(tol = 1E-20), cW = 0.1){
  inWW <- rep(TRUE, ncol(sample))
  ldenfun <- function(sample, theta){ #here theta is km
    return(dmovMF(sample, theta, log = TRUE))
  }
  if (method == "smfull"){
    if (is.null(km)){fixedtheta <- !is.na(km)}
    else {fixedtheta <- rep(FALSE, ncol(sample))}
    starttheta <- vMF_full(sample, km = km, control = control)
    estimator <- function(Y, weights, starttheta, fixedtheta){
      out <- vMF_full(sample, km = starttheta * fixedtheta,
                      control = control, w=w, starttheta)
      return(out$km)
    }
  } else if (method == "Mardia"){
    stopifnot(is.null(km))
    fixedtheta <- rep(FALSE, ncol(sample))
    starttheta <- vMF_Mardia(sample, control = control)
    estimator <- function(Y, weights, starttheta, fixedtheta){
      out <- vMF_Mardia(sample, control = control, w=w)
      return(out$km)
    }
  }

  est <- windham_raw(prop = m$sample,
                     cW = cW,
                     ldenfun = ldenfun,
                     estimatorfun = estimator,
                     starttheta = km,
                     fixedtheta = fixedtheta,
                     inWW = inWW,
                     originalcorrectionmethod = TRUE)
  return(est)
}

vMF_Mardia <- function(sample, control = list(tol = 1E-20), w = rep(1, nrow(sample))){
  mu <- apply(sample, MARGIN = 2, weighted.mean, w)
  mu <- mu/sqrt(sum(mu^2))
  Rtrans <- Directional::rotation(mu, c(1, rep(0, length(mu) - 1)))
  samplestd <- sample %*% t(Rtrans)
  # check: mustd <- colMeans(samplestd); mustd <- mustd / sqrt(sum(mustd^2))

  # do estimate, where all but the first component of theta are fixed at zero
  # because kappa * e1 = (kappa, 0, 0, 0, ...)
  p <- ncol(sample)
  intheta <- c(NA, rep(0, p - 1))
  tapes <- buildsmotape("Snative", "vMF",
                        rep(1, p)/sqrt(p), intheta,
                        weightname = "ones",
                        verbose = FALSE)
  out <- smest(tapes$smotape, rep(0.1, sum(is.na(intheta))), samplestd, control = control, w = w)
  return(list(
    k = out$par,
    m = mu,
    km = out$par * mu,
    SE = list(k = out$SE),
    sminfo = out
  ))
}

vMF_full <- function(sample, km = NULL, control = list(tol = 1E-20), w = w, starttheta = NULL){
  p <- ncol(sample)
  if (is.null(km)){
    km <- rep(NA, p)
  }
  stopifnot(length(km) == p)

  intheta <- km
  tapes <- buildsmotape("Snative", "vMF",
                        rep(1, p)/sqrt(p), intheta,
                        weightname = "ones",
                        verbose = FALSE)
  if (is.null(starttheta)){starttheta <- rep(0.1, sum(is.na(intheta)))}
  else {starttheta <- starttheta[is.na(intheta)]}
  out <- smest(tapes$smotape, starttheta, sample, control = control, w=w)
  theta <- intheta
  theta[is.na(intheta)] <- out$par

  SE <- intheta * 0
  SE[is.na(intheta)] <- out$SE
  return(list(
    km = theta,
    k = sqrt(sum(theta^2)),
    m = theta / sqrt(sum(theta^2)),
    SE = list(km = SE),
    sminfo = out
  ))
}
