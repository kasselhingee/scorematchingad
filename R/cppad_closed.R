#' @title Score Matching Estimator for Quadratic-Form Score-Matching Objectives
#' @description 
#' Returns the vector when the gradient of the score matching objective is zero for any quadratic-form tape.
#' Also estimates standard errors and covariance.
#' @param Yapproxcentres A matrix of Taylor approximation centres for rows of Y that require approximation. `NA` for rows that do not require approximation.
#' @param smotape A tape of a score matching objective function that has *quadratic form*. Test for quadratic form using [`testquadratictape()`].
#' The `smotape`'s independent variables are assumed to be the model parameters to fit
#' and the `smotape`'s dynamic parameter is a (multivariate) measurement.
#' @param Y A matrix of multivariate observations. Each row is an observation.
#' @param w Weights for each observation.
#' @param approxorder The order of Taylor approximation to use.
#' @details
#' Standard errors are only computed when the weights are constant, and use the Godambe information matrix (aka sandwich method) [*would like a good reference for this here*].
#' The sensitivity matrix \eqn{G} is estimated as
#' the negative of the average over the Hessian of `smotape` evaluated at each observation in `Y`.
# \deqn{\hat{G(\theta)} = \hat{E} -H(smo(\theta;Y))),}
# where \eqn{smo} is the score matching objective function represented by `smotape`,
# \eqn{H} is the Hessian with respect to \eqn{\theta}, which is constant for quadratic-form functions,
# 
#' The variability matrix \eqn{J} is then estimated as
#' the sample covariance (denominator of \eqn{n-1}) of the gradiant of `smotape` evaluated at each of the observations in `Y` for the estimated \eqn{\theta}.
# \deqn{\hat{J}(\theta) = var(grad(w smo(\theta;Y))),}

#' The variance of the estimator is then estimated as
#' \eqn{G^{-1}JG^{-1}/n,}
# \deqn{\hat{G}(\theta)^{-1}\hat{J}(\theta)\hat{G}(\theta)^{-1}/n,}
#' where `n` is the number of observations.

cppad_closed <- function(smotape, Y, Yapproxcentres = NA * Y, 
                         w = rep(1, nrow(Y)),
                         approxorder = 10){
  stopifnot(nrow(Y) == length(w))
  stopifnot(nrow(Y) == nrow(Yapproxcentres))
  
  parts <- quadratictape_parts(smotape, tmat = Y,
                               tcentres = Yapproxcentres,
                               approxorder = approxorder)

  # weight parts
  parts$offset <- parts$offset * w
  parts$Hessian <- parts$Hessian * w
  offset <- colSums(parts$offset)
  Hess <- colSums(parts$Hessian)
  Hess <- matrix(Hess, ncol = sqrt(ncol(parts$Hessian)))
  invHess <- solve(Hess)
  root <- drop(-1 * invHess %*% offset)

  # compute SEs
  SE <- NA
  covar <- NA
  if (all(w[[1]] == w)){
    covar <- cppad_closed_estvar(Y, root, parts$offset, parts$Hess)
    SE <- sqrt(diag(estcovar))
  }
  return(list(
    est = root,
    Hessian = Hess,
    offset = offset,
    SE = SE,
    covar = covar
  ))
}

# function for getting standard errors for the above function
cppad_closed_estvar <- function(Y, theta, offsets, Hesss){
  sens <- - matrix(colMeans(Hesss), ncol = sqrt(ncol(Hesss)))
  sensinv <- solve(sens)

  # compute gradients
  grads <- lapply(1:nrow(offsets), function(i){
    drop(matrix(Hesss[i, ], ncol = ncol(offsets)) %*% theta + 
      t(offsets[i, , drop = FALSE]))
  })
  grads <- do.call(rbind, grads)
  variability <- cov(grads)

  Ginfinv <- sensinv %*% variability %*% sensinv / nrow(offsets) #inverse of the Godambe information matrix, also called the sandwich information matrix
  return(Ginfinv) 
}


