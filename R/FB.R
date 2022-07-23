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
#' @param km Optional. A vector of same length as the dimension, representing the elements of the Fisher component (concentraction * vector).
#' If supplied, the non-NA elements are fixed.
#' @param A Optional. The Bingham matrix. If supplied the non-NA elements of the Bingham matrix is fixed.
#' The final diagonal of `A` here must be NA as the software calculates this value to ensure the trace is zero.
#' @details
#' Warning: the score matching estimate appears to converge very slowly for the Fisher-Bingham distribution.
#' Even with a million simulated measurements,
#'  the gradient of the score matching objective at the true theta is
#'  is of size 0.001, which is substantially non-zero.
#' @export
FB <- function(sample, km = NULL, A = NULL, control = default_Rcgmin()){
  p <- ncol(sample)
  if (is.null(A)){
    A <- matrix(NA, nrow = p, ncol = p)
  }
  stopifnot(all(dim(A) == c(p, p)))
  if (!is.na(A[p,p])){stop("The final diagonal element of matrix A cannot be fixed in this software. Please consider reordering your dimensions.")}
  if (is.null(km)){
    km <- rep(NA, p)
  }

  intheta <- FB_mats2theta(1, km, A)

  tapes <- buildsmotape("Snative", "FB",
               rep(1, p)/sqrt(p), intheta,
               weightname = "ones",
               verbose = FALSE)

  sminfo <- cppadest(tapes$smotape, rep(0.1, sum(is.na(intheta))), sample,
               control = control)
  theta <- intheta
  theta[is.na(intheta)] <- sminfo$par
  thetamat <- FB_theta2mats(theta)
  SE <- intheta * 0
  SE[is.na(intheta)] <- sminfo$SE
  SE <- FB_theta2mats(SE, isSE = TRUE)
  return(c(thetamat, list(sminfo = sminfo, SE = SE)))
}

# non-normalised density function
qdFB <- function(u, k, m, A){
  exp(u %*% A %*% u + k * m %*% u)
}

FB_mats2theta <- function(k, m, A){
  Binghamtheta <- Bingham_Amat2theta(A)
  theta <- c(Binghamtheta, k * m)
  return(theta)
}

FB_theta2mats <- function(theta, isSE = FALSE){
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
  km <- theta[seq.int(p-1 + (p - 1) * p/2 + 1, length.out = p)]
  if (!isSE){
    A[p,p] <- -sum(diag(A)[1:(p-1)]) #trace is 0 constraint
    k <- sqrt(sum(km^2))
    m <- km/k
    return(list(
      k = k,
      m = m,
      km = km,
      A = A
    ))
  } else {
    A[p,p] <- NA #trace is 0 constraint
    return(list(
      km = km,
      A = A
    ))
  }
}
