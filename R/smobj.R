#' @title Compute value of score matching objective for multiple observations
#' @description The value of the score matching objective for a given beta
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
#' @export
smobj <- function(smofun, theta, utabl,
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
  if (is.null(w)){scmo <- mean(unlist(sc_perpt))}
  else {scmo <- weighted.mean(unlist(sc_perpt), w=w)}
  if (stopifnan && is.nan(scmo)){stop("smobj() generates a NaN")}
  return(scmo)
}

smobj_b <- function(theta, ...){
  args = list(...)
  args2 <- args[names(args) %in% formalArgs(smobj)]
  do.call(smobj, c(list(theta = theta), args2))
}

#' @describeIn smobj The gradient of the score matching objective function at given beta
#' @export
smobjgrad <- function(smofun, theta, utabl,
                      Jsmofun_u = NULL, uboundary = NULL, boundaryapprox = NULL, approxorder = NULL,
                      stopifnan = FALSE, w = NULL){
  grad_perpt <- smobjgrad_perpt(smofun, theta, utabl,
                  Jsmofun_u = Jsmofun_u, uboundary = uboundary, boundaryapprox = boundaryapprox,
                  approxorder = approxorder)
  if (is.null(w)){grad <- colMeans(do.call(rbind, grad_perpt))}
  else {grad <- apply(do.call(rbind, grad_perpt), MARGIN = 2, weighted.mean, w=w)}
  if (stopifnan && any(is.nan(grad))){stop("smobjgrad() generates a NaN")}
  return(grad)
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

smobjgrad_b <- function(theta, ...){
  args = list(...)
  args2 <- args[names(args) %in% formalArgs(smobjgrad)]
  do.call(smobjgrad, c(list(theta = theta), args2))
}

#' @describeIn smobj Hessian of the score matching objective
#' @export
smobjhess <- function(smofun, theta, utabl,
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

  if (is.null(w)){hess <- colMeans(do.call(rbind, hess_perpt))}
  else {hess <- apply(do.call(rbind, hess_perpt), MARGIN = 2, weighted.mean, w=w)}
  dim(hess) <- rep(length(theta), 2)
  if (stopifnan && any(is.nan(hess))){stop("smobjhess() generates a NaN")}
  return(hess)
}
