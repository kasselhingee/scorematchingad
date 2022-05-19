#' @title Score matching estimates for the Bingham distribution
#' @param sample Samples from the Bingham distribution in cartesian Rd coordinates. Each row is a measurement.
#' @param control Control parameters passed to `Rcgmin::Rcgmin()`

#' @description The Bingham distribution as described by Mardia et al 2016.
#'  Book Mardia and Jupp 2000 would likely be a better reference.
#' @examples
#' p <- 4
#' A <- matrix(NA, ncol = p, nrow = p)
#' A[upper.tri(A)] <- runif(sum(upper.tri(A)))
#' A[lower.tri(A)] <- t(A)[lower.tri(A)]
#' diag(A) <- c(runif(p-1), NA)
#' A[p,p] <- -sum(diag(A)[1:(p-1)]) #trace is 0 constraint
#' sample <- rBingham(100, A)
#' @export


#' @describeIn Bingham Uses `Directional::rbingham` to simulate from the Bingham distribution.
#' @param A is the parameter matrix for the Bingham distribution. The log density is proportional to `t(u) * A * u`
#' @param n The number of samples to generate.
#' @export
rBingham <- function(n, A){
  if (is.matrix(A)){
    if (dim(A)[[1]] == dim(A)[[2]]){
      sample <- Directional::rbingham(n, -A) # Fallaize and Kypraios 2018 use a '-' in the exponent for the Bingham distribution
    }
  } else {
    A <- Bingham_theta2Amat(A)
    sample <- Directional::rbingham(n, -A)
  }
  return(sample)
}

Bingham_full <- function(sample, control = list(tol = 1E-20)){
  pman <- pmanifold("Snative")
  p <- ncol(sample)
  thetatape <- rep(0.1, p-1 + (p - 1) * p/2)
  lltape <- ptapell(rep(1, p) / sqrt(p), thetatape,
                    llname = "Bingham", pman,
                    fixedtheta = rep(FALSE, length(thetatape)), verbose = FALSE)
  smotape<- ptapesmo(rep(1, p) / sqrt(p), thetatape,
                     pll = lltape, pman = pman, "ones", acut = 1, verbose = FALSE) #tape of the score function
  out <- smest(smotape, thetatape, sample, control = control)
  SE = Bingham_theta2Amat(out$SE)
  SE[p, p] <- NA
  return(list(
    A = Bingham_theta2Amat(out$par),
    SE = SE,
    sminfo = out
  ))
}

Bingham_Amat2theta <- function(A){
  p <- ncol(A)
  stopifnot(isSymmetric(A))
  stopifnot(abs(sum(diag(A))) < 1E-8)
  theta <- c(diag(A)[1:(p-1)], A[upper.tri(A)])
  return(theta)
}

Bingham_theta2Amat <- function(theta){
  #length of theta is p-1 + (p - 1) * p/2 = p2/2 - p/2 + p - 1 = p2/2 + p/2 - 1
  # => 0 = p2/2 + p/2 - 1 - length
  # => 0 = p2 + p - 2(1+length)
  # => p = (-1 +/- sqrt(1 + 4 * 2 * (1+length))) / 2
  p = (-1 + sqrt(1 + 4 * 2 * (1+length(theta)))) / 2  #the '-' alternative is always negative because sqrt(1 +epsilon) > 1
  A <- matrix(NA, ncol = p, nrow = p)
  diag(A) <- c(theta[1:(p-1)], NA)
  A[upper.tri(A)] <- theta[p:length(theta)]
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #trace is 0 constraint
  return(A)
}

