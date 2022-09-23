#' @noRd
#' @title Compute value of score matching objective for multiple observations
#' @description Computes the score matching objective given theta
#' @param smofun A tape of the score matching objective calculation
#' @param theta A parameter set
#' @param utabl A matrix of observations, each row being an observation.
#' @param uboundary Optional. A matrix of observations close to or on the boundary.
#' The contribution of the score matching objective (and corresponding derivatives) of these measurements
#' will be approximated using Taylor approximation from the corresponding `boundaryapprox` matrix.
#' @param boundaryapprox Optional. A matrix of locations in the manifold associated with each measurement in `boundary`.
#' For best results these locations should be further from the manifold boundary and close to their corresponding measurements.
#' Taylor approximation around the rows of `boundaryapprox` will be used to approximate the score matching objective for these measurements.
#' @param stopifnan If TRUE function will create an error when the return is not a number.
#' @return The functions `*_sum` all return a value with an attribute 'normaliser' that is the normalising constant needed to turn the result into an average. The functions `smobj()`, `smobjgrad()`, `smobjhess()` are wrapper of the `*_sum` functions for easier testing and UI.
# @export
smobj <- function(...){ #smobj average of observations (i.e. smobj_sum / sum(w))
  args <- list(...)
  asum <- do.call(smobj_sum, args)
  avg <- asum / attr(asum, "normaliser")
  attr(avg, "normaliser") <- NULL
  return(avg)
}

# @describeIn smobj The sum of score-matching objective value (no division by n - useful for optimisation where floating point issues can accumulate)
smobj_sum <- function(smofun, theta, utabl,
                  smofun_u = NULL, uboundary = NULL, boundaryapprox = NULL, approxorder = NULL,
                  stopifnan = FALSE, w = NULL){
  sc_perpt_interior <- c()
  if (nrow(utabl) > 0){
    sc_perpt_interior <- lapply(1:nrow(utabl), function(i){
      scobj <- pForward0(smofun, theta, utabl[i,])
      return(scobj)
    })
  }

  sc_perpt_boundary <- c()
  if (any(isTRUE(nrow(uboundary) > 0), isTRUE(nrow(boundaryapprox) > 0))){
    stopifnot(all(!is.null(smofun_u), isTRUE(nrow(uboundary) > 0), isTRUE(nrow(boundaryapprox) > 0), !is.null(approxorder)))
    stopifnot(nrow(uboundary) == nrow(boundaryapprox))
    sc_perpt_boundary <- lapply(1:nrow(uboundary), function(i){
      scobj <- pTaylorApprox(smofun_u, uboundary[i,], boundaryapprox[i, ], theta, order = approxorder)
      return(scobj)
    })
  }
  sc_perpt <- c(sc_perpt_interior, sc_perpt_boundary)
  if (is.null(w)){
     nscmo <- sum(unlist(sc_perpt))
     attr(nscmo, "normaliser") <- length(sc_perpt)
  } else {
     stopifnot(length(w) == length(sc_perpt))
     nscmo <- sum(unlist(sc_perpt) * w)
     attr(nscmo, "normaliser") <- sum(w)}
  if (stopifnan && is.nan(nscmo)){stop("smobj() generates a NaN")}
  return(nscmo)
}

smobj_sum_b <- function(theta, ...){
  args = list(...)
  args2 <- args[names(args) %in% formalArgs(smobj_sum)]
  do.call(smobj_sum, c(list(theta = theta), args2))
}

# @describeIn smobj The gradient of the score matching objective function at given beta
# @export
smobjgrad <- function(...){ #smobjgrad average of observations (i.e. smobjgrad_sum / sum(w))
  args <- list(...)
  asum <- do.call(smobjgrad_sum, args)
  avg <- asum / attr(asum, "normaliser")
  attr(avg, "normaliser") <- NULL
  return(avg)
}

smobjgrad_sum <- function(smofun, theta, utabl,
                      Jsmofun_u = NULL, uboundary = NULL, boundaryapprox = NULL, approxorder = NULL,
                      stopifnan = FALSE, w = NULL){
  grad_perpt <- smobjgrad_perpt(smofun, theta, utabl,
                  Jsmofun_u = Jsmofun_u, uboundary = uboundary, boundaryapprox = boundaryapprox,
                  approxorder = approxorder)
  grad_perpt <- do.call(rbind, grad_perpt)
  if (is.null(w)){
    ngrad <- colSums(grad_perpt)
    attr(ngrad, "normaliser") <- nrow(grad_perpt)
  } else {
    stopifnot(length(w) == nrow(grad_perpt))
    ngrad <- colSums(grad_perpt * w)
    attr(ngrad, "normaliser") <- sum(w)
  }
  if (stopifnan && any(is.nan(ngrad))){stop("smobjgrad() generates a NaN")}
  return(ngrad)
}

smobjgrad_perpt <- function(smofun, theta, utabl,
                      Jsmofun_u = NULL, uboundary = NULL, boundaryapprox = NULL, approxorder = NULL){
  grad_perpt_interior <- list()
  if (nrow(utabl) > 0){
    grad_perpt_interior <- lapply(1:nrow(utabl), function(i){
      grad <- pJacobian(smofun, theta, utabl[i,])
      return(grad)
    })
  }

  grad_perpt_boundary <- list()
  if (any(isTRUE(nrow(uboundary) > 0), isTRUE(nrow(boundaryapprox) > 0))){
    stopifnot(all(!is.null(Jsmofun_u), isTRUE(nrow(uboundary) > 0), isTRUE(nrow(boundaryapprox) > 0), !is.null(approxorder)))
    stopifnot(nrow(uboundary) == nrow(boundaryapprox))
    grad_perpt_boundary <- lapply(1:nrow(uboundary), function(i){
      grad <- pTaylorApprox(Jsmofun_u, uboundary[i,], boundaryapprox[i, ], theta, order = approxorder)
      return(grad)
    })
  }
  grad_perpt <- c(grad_perpt_interior, grad_perpt_boundary)
  return(grad_perpt)
}

smobjgrad_sum_b <- function(theta, ...){
  args = list(...)
  args2 <- args[names(args) %in% formalArgs(smobjgrad_sum)]
  do.call(smobjgrad_sum, c(list(theta = theta), args2))
}

# @describeIn smobj Hessian of the score matching objective
# @export
smobjhess <- function(...){ #smobjgrad average of observations (i.e. smobjgrad_sum / sum(w))
  args <- list(...)
  asum <- do.call(smobjhess_sum, args)
  avg <- asum / attr(asum, "normaliser")
  attr(avg, "normaliser") <- NULL
  return(avg)
}
smobjhess_sum <- function(smofun, theta, utabl,
                      Hsmofun_u = NULL, uboundary = NULL, boundaryapprox = NULL, approxorder = NULL,
                      stopifnan = FALSE, w = NULL){
  hess_perpt_interior <- list()
  if (nrow(utabl) > 0){
    hess_perpt_interior <- lapply(1:nrow(utabl), function(i){
      hess <- pHessian(smofun, theta, utabl[i,])
      return(hess)
    })
  }

  hess_perpt_boundary <- list()
  if (any(isTRUE(nrow(uboundary) > 0), isTRUE(nrow(boundaryapprox) > 0))){
    stopifnot(all(!is.null(Hsmofun_u), isTRUE(nrow(uboundary) > 0), isTRUE(nrow(boundaryapprox) > 0), !is.null(approxorder)))
    stopifnot(nrow(uboundary) == nrow(boundaryapprox))
    hess_perpt_boundary <- lapply(1:nrow(uboundary), function(i){
      hess <- pTaylorApprox(Hsmofun_u, uboundary[i,], boundaryapprox[i, ], theta, order = approxorder)
      return(hess)
    })
  }

  hess_perpt <- c(hess_perpt_interior, hess_perpt_boundary)
  hess_perpt <- do.call(rbind, hess_perpt)
  if (is.null(w)){
    nhess <- colSums(hess_perpt)
    attr(nhess, "normaliser") <- nrow(hess_perpt)
  } else {
    stopifnot(length(w) == nrow(hess_perpt))
    nhess <- colSums(hess_perpt * w)
    attr(nhess, "normaliser") <- sum(w)}
  dim(nhess) <- rep(length(theta), 2)
  if (stopifnan && any(is.nan(nhess))){stop("smobjhess() generates a NaN")}
  return(nhess)
}

