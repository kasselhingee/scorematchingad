#' @title Score matching estimates for the Bingham distribution
#' @param sample Samples from the Bingham distribution in cartesian Rd coordinates. Each row is a measurement.
#' @param control Control parameters passed to `Rcgmin::Rcgmin()`
#' @param method The estimating method, either `smfull` for score matching estimates for all parameters
#'  or `Mardia` for the Mardia et al hybrid method. `hybrid` will also select the Mardia et al method.

#' @description The Bingham distribution as described by Mardia et al 2016.
#'  Book Mardia and Jupp 2000 would likely be a better reference.
#' @examples
#' p <- 4
#' A <- rsymmetricmatrix(p)
#' A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
#' sample <- rBingham(100, A)
#'
#' Bingham(sample, method = "Mardia")
#' @export
Bingham <- function(sample, method = "smfull", control = list(tol = 1E-20)){
  if (method == "smfull"){out = Bingham_full(sample, control = control)}
  if (method == "Mardia"){out <- Bingham_Mardia(sample, control = control)}
  if (method == "hybrid"){out <- Bingham_Mardia(sample, control = control)}
  return(out)
}

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
  A = Bingham_theta2Amat(out$par)
  SE = Bingham_theta2Amat(out$SE)
  SE[p, p] <- NA

  A_es <- eigen(A)
  return(list(
    A = A,
    A_SE = SE,
    Gamma = A_es$vectors,
    Lambda = A_es$values,
    sminfo = out
  ))
}

Bingham_Mardia <- function(sample, control = list(tol = 1E-20)){
  Tmat <- 1/nrow(sample) * t(sample) %*% sample
  Tmat_es <- eigen(Tmat)
  Gammahat <- Tmat_es$vectors
  samplestd <- sample %*% Gammahat

  pman <- pmanifold("Snative")
  p <- ncol(samplestd)
  thetatape <- c(rep(0.1, p-1), rep(0, (p - 1) * p/2)) # A for standardised data is a diagonal matrix
  lltape <- ptapell(rep(1, p) / sqrt(p), thetatape,
                    llname = "Bingham", pman,
                    fixedtheta = c(rep(FALSE, p-1), rep(TRUE,  (p - 1) * p/2)), verbose = FALSE)
  smotape<- ptapesmo(rep(1, p) / sqrt(p), thetatape[1:(p-1)],
                     pll = lltape, pman = pman, "ones", acut = 1, verbose = FALSE) #tape of the score function

  sm <- smest(smotape, thetatape[1:(p-1)], samplestd, control = control)
  Lambda <- c(sm$par, -sum(sm$par))
  return(list(
    Lambda = Lambda,
    Lambda_SE = c(sm$SE, NA),
    Gamma = Gammahat,
    A = Gammahat %*% diag(Lambda) %*% t(Gammahat),
    sminfo = sm
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

