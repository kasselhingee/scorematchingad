#' @title Score matching estimates for the Bingham distribution
#' @param sample A matrix of observations in Cartesian coordinates. Each row is a (multivariate) measurement.
#' @param control Control parameters passed to [`Rcgmin::Rcgmin()`].
#' @param A For full score matching only: if supplied, then NA elements of `A` are estimated and the other elements are fixed.
#' @param method The estimating method, either "smfull" for score matching estimates for all parameters
#'  or "Mardia" for the \insertCite{mardia2016sc}{scorecompdir} hybrid estimator.


#' @description
#' Score matching estimator of the Bingham matrix.
#' @details
#' The Bingham distribution , which has a density proportional to
#' \deqn{\exp(z^T A z),}
#' where \eqn{A} is a symmetric matrix and the trace (sum of the diagonals) of \eqn{A} zero for identifiability \insertCite{@p181, @mardia2000di}{scorecompdir}.
#'
#' The function `Bingham()` calculates the final diagonal element of \eqn{A} from the sum of the other diagonal elements to ensure that the trace of \eqn{A} is zero. All other elements of \eqn{A} are estimated using score matching.
#'
#' The estimating method by \insertCite{mardia2016sc;textual}{scorecompdir} first calculates the maximum-likelihood estimate of the eigenvectors \eqn{G} of \eqn{A}. 
#' The observations `Y` are then standardised to `Y`\eqn{G}. 
#' After this standardisation, the non-diagonal elements of \eqn{A} are zero.
#' The diagonal elements are estimated using score matching.
#' See \insertCite{mardia2016sc}{scorecompdir} for details.
#' @examples
#' p <- 4
#' A <- rsymmetricmatrix(p)
#' A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
#' sample <- rBingham(100, A)
#'
#' Bingham(sample, method = "Mardia")
#' @family Mardia hybrid estimators
#' @export
Bingham <- function(sample, A = NULL, method = "smfull", control = default_Rcgmin()){
  if (method == "smfull"){
    out <- Bingham_full(sample, A = A, control = control)}
  if (method %in% c("Mardia", "hybrid")){
    stopifnot(is.null(A))
    out <- Bingham_Mardia(sample, control = control)
    }
  return(out)
}

#' @noRd
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

Bingham_full <- function(sample,  A = NULL, control = default_Rcgmin()){
  p <- ncol(sample)
  if (is.null(A)){
    A <- matrix(NA, nrow = p, ncol = p)
  }
  stopifnot(all(dim(A) == c(p, p)))
  if (!is.na(A[p,p])){stop("The final diagonal element of matrix A cannot be fixed in this software. Please consider reordering your dimensions.")}

  intheta <- Bingham_Amat2theta(A)

  utape <- rep(1, p) / sqrt(p)
  tapes <- buildsmotape("Snative", "Bingham",
                           utape, intheta,
                           weightname = "ones")
  out <- cppadest(tapes$smotape, rep(0.1, sum(is.na(intheta))), sample, control = control)
  theta <- intheta
  theta[is.na(intheta)] <- out$par
  tSE <- intheta * 0
  tSE[is.na(intheta)] <- out$SE
  A = Bingham_theta2Amat(theta)
  SE = Bingham_theta2Amat(tSE)
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

Bingham_Mardia <- function(sample,  control = default_Rcgmin()){
  Tmat <- 1/nrow(sample) * t(sample) %*% sample
  Tmat_es <- eigen(Tmat)
  Gammahat <- Tmat_es$vectors
  samplestd <- sample %*% Gammahat

  p <- ncol(samplestd)
  A <- matrix(0, nrow = p, ncol = p)
  diag(A) <- NA
  stopifnot(all(dim(A) == p))
  intheta <- Bingham_Amat2theta(A)
  tapes <- buildsmotape("Snative", "Bingham",
                        rep(1, p) / sqrt(p), intheta,
                        weightname = "ones")

  sm <- cppadest(tapes$smotape, seq(length.out = sum(is.na(intheta))), samplestd, control = control)
  theta <- intheta
  theta[is.na(intheta)] <- sm$par
  Astd <- Bingham_theta2Amat(theta)
  Lambda <- diag(Astd)

  SE <- intheta * 0
  SE[is.na(intheta)] <- sm$SE
  SE <- Bingham_theta2Amat(SE)
  SE[p, p] <- NA#final NA because final diagonal element is not estimated directly
  return(list(
    Lambda = Lambda,
    Lambda_SE = diag(SE),
    Gamma = Gammahat,
    A = Gammahat %*% diag(Lambda) %*% t(Gammahat),
    sminfo = sm
  ))
}

Bingham_Amat2theta <- function(A){
  p <- ncol(A)
  stopifnot(isSymmetric(A))
  if(isTRUE(abs(sum(diag(A))) > 1E-8)){warning("Trace of A is not zero, final diagonal element of A will be altered.")}
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

