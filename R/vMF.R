#' @title Score matching estimate for the von-Mises Fisher distribution
#' @param sample Samples from the vMF distribution in cartesian Rd coordinates. Each row is a measurement.
#' @param control Control parameters passed to `Rcgmin::Rcgmin()`
#' @examples
#' sample <- movMF::rmovMF(100, 3 * c(1, -1) / sqrt(2))
#' vMF(sample, "smfull")
#' vMF(sample, "Mardia")
#' @export
vMF <- function(sample, method = "smfull", control = list(tol = 1E-20)){
  out <- NULL
  if (method == "smfull"){out <- vMF_full(sample, control = list(tol = 1E-20))}
  if (method == "Mardia"){
    out <- vMF_Mardia(sample, control = list(tol = 1E-20))
  }
  return(out)
}

vMF_Mardia <- function(sample, control = list(tol = 1E-20)){
  mu <- colMeans(sample)
  mu <- mu/sqrt(sum(mu^2))
  R <- Directional::rotation(mu, c(1, rep(0, length(mu) - 1)))
  samplestd <- t(R %*% t(sample))
  # check: mustd <- colMeans(samplestd); mustd <- mustd / sqrt(sum(mustd^2))
  # do estimate, where all but the first component of theta are fixed at zero
  # because kappa * e1 = (kappa, 0, 0, 0, ...)
  pman <- pmanifold("Snative")
  p <- ncol(sample)
  lltape <- ptapell(rep(1, p) / p, c(1, rep(0, p - 1)), llname = "vMF", pman,
                    fixedtheta = c(FALSE, rep(TRUE, p - 1)), verbose = FALSE)
  smotape<- ptapesmo(rep(1, p) / p, 1,
                     pll = lltape, pman = pman, "ones", acut = 1, verbose = FALSE) #tape of the score function
  out <- smest(smotape, 1, sample, control = control)
  return(list(
    kappa = out$par,
    mu = mu,
    sminfo = out
  ))
}

vMF_full <- function(sample, control = list(tol = 1E-20)){
  pman <- pmanifold("Snative")
  p <- ncol(sample)
  lltape <- ptapell(rep(1, p) / p, rep(0.1, p), llname = "vMF", pman,
                    fixedtheta = rep(FALSE, p), verbose = FALSE)
  smotape<- ptapesmo(rep(1, p) / p, rep(0.1, p),
                     pll = lltape, pman = pman, "ones", acut = 1, verbose = FALSE) #tape of the score function
  out <- smest(smotape, rep(0.1, p), sample, control = control)
  return(list(
    kappa = sqrt(sum(out$par^2)),
    mu = out$par / sqrt(sum(out$par^2)),
    sminfo = out
  ))
}
