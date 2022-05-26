#' @title Estimate Fisher-Bingham distribution
#' @param control A list of control parameters passed to `Rcgmin::Rcgmin()`.
#' @examples
#' p <- 3
#' A <- rsymmetricmatrix(p, -10, 10)
#' A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
#' m <- runif(p, -10, 10)
#' m <- m / sqrt(sum(m^2))
#' sample <- rFB(1000, 2, m, A)
#' FB(sample)
#' @param sample An array of samples, each row an individual sample.
#' @export
FB <- function(sample, control = list(tol = 1E-20)){
  p <- ncol(sample)
  ltheta <- p-1 + (p - 1) * p/2 + p
  thetatape <- seq.int(1, length.out = ltheta)
  pman <- pmanifold("Snative")

  lltape <- ptapell(sample[1,], thetatape, llname = "FB", pman,
                    fixedtheta = rep(FALSE, ltheta), verbose = FALSE)
  smotape <- ptapesmo(sample[1,], thetatape,
                      lltape, pman, "ones", 1, verbose = FALSE)
  sminfo <- smest(smotape, thetatape, sample,
               control = control)
  thetamat <- FB_theta2mats(sminfo$par)
  return(c(thetamat, list(sminfo = sminfo)))
}

# non-normalised density function
qdFB <- function(u, k, m, A){
  exp(u %*% A %*% u + k * m %*% u)
}

FB_mats2theta <- function(k, m, A){
  p <- ncol(A)
  stopifnot(isSymmetric(A))
  stopifnot(abs(sum(diag(A))) < 1E-8)
  theta <- c(diag(A)[1:(p-1)], A[upper.tri(A)], k * m)
  return(theta)
}

FB_theta2mats <- function(theta){
  #length of theta is:
  # l = p-1 + (p - 1) * p/2 + p
  # = p^2/2 + p * (1 - 1/2 + 1) - 1
  # = p^2/2 + p * (1.5) - 1
  # so
  # 0 = p^2 + 3*p - (2 + 2*l)
  # => p = (-3 pm sqrt(3^2 + 4(2+2l)))/2
  #  = (-3 pm sqrt(9 + 8 + 8l)))/2
  #  = (-3 pm sqrt(17 + 8l)))/2
  # = (-3 + sqrt(17 + 8l)))/2  # because sqrt(17 + 8l) >= 3
  p <- (-3 + sqrt(17 + 8 * length(theta)))/2
  A <- matrix(NA, ncol = p, nrow = p)
  diag(A) <- c(theta[1:(p-1)], NA)
  A[upper.tri(A)] <- theta[seq.int(p, length.out = (p - 1) * p/2)]
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #trace is 0 constraint
  m <- theta[seq.int(p-1 + (p - 1) * p/2 + 1, length.out = p)]
  k <- sqrt(sum(m^2))
  m <- m/k
  return(list(
    k = k,
    m = m,
    A = A
  ))
}
