#' @title Score matching estimate for the von-Mises Fisher distribution
#' @param sample Samples from the vMF distribution in cartesian Rd coordinates. Each row is a measurement.
#' @param control Control parameters passed to `Rcgmin::Rcgmin()`
#' @examples
#' sample <- movMF::rmovMF(100, 3 * c(1, -1) / sqrt(2))
#' vMF(sample)
#' @export
vMF <- function(sample, control = list(tol = 1E-20)){
  pman <- pmanifold("Snative")
  p <- ncol(sample)
  lltape <- ptapell(rep(1, p) / p, rep(0.1, p), llname = "vMF", pman,
                    fixedtheta = rep(FALSE, p), verbose = FALSE)
  smotape<- ptapesmo(rep(1, p) / p, rep(0.1, p),
                     pll = lltape, pman = pman, "ones", acut = 1, verbose = FALSE) #tape of the score function
  out <- smest(smotape, rep(0.1, p), sample, control = list(tol = 1E-10))
  return(out)
}
