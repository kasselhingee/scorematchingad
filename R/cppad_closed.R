#' @title Score Matching Estimator for Quadratic-Form Score-Matching Objective Functions
#' @description 
#' For a `CppAD` tape of a quadratic-form score matching objective function, calculates the vector of parameters such that the gradient of the score matching objective is zero.
#' Also estimates standard errors and covariance.
#' Many score matching estimators have an objective function that has a quadratic form.
#' @param Yapproxcentres A matrix of Taylor approximation centres for rows of Y that require approximation. `NA` for rows that do not require approximation.
#' @param smotape A `CppAD` tape of a score matching objective function that has *quadratic form*. Test for quadratic form using [`testquadratictape()`].
#' The `smotape`'s independent variables are assumed to be the model parameters to fit
#' and the `smotape`'s dynamic parameter is a (multivariate) measurement.
#' @param Y A matrix of multivariate observations. Each row is an observation.
#' @param w Weights for each observation.
#' @param approxorder The order of Taylor approximation to use.
#' @details
#' When the score matching function is of quadratic form, then the gradient is zero at \eqn{H^{-1}b}{solve(H) %*% b},
#' where \eqn{H} is the average of the Hessian of the score matching function evaluated at each measurement and 
#' \eqn{b} is the average of the gradient offset (see [`quadratictape_parts()`]) evaluated at each measurement. 
#' Both the Hessian and the gradient offset are constant with respect to the model parameters for quadratic-form score matching equations.
#'
#' Standard errors are only computed when the weights are constant, and use the Godambe information matrix (aka sandwich method) (*would like a good reference for this here*).
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
#' @family {General score matching}
#' @export
cppad_closed <- function(smotape, Y, Yapproxcentres = NA * Y, 
                         w = rep(1, nrow(Y)),
                         approxorder = 10){
  stopifnot(nrow(Y) == length(w))
  stopifnot(nrow(Y) == nrow(Yapproxcentres))
  
  parts <- quadratictape_parts(smotape, tmat = Y,
                               tcentres = Yapproxcentres,
                               approxorder = approxorder)

  # weight parts
  partsT <- lapply(parts, wcolSums, w = w)
  offset <- partsT$offset
  Hess <- partsT$Hessian
  Hess <- matrix(Hess, ncol = sqrt(ncol(parts$Hessian)))
  invHess <- tryCatch(solve(Hess),
               error = function(e) {
     if (grepl("system.*singular", e)){
          stop(paste("Hessian of the score-matching objective function is not invertible.", e))
        } else {stop(e)}
     })
  root <- drop(-1 * invHess %*% offset)

  # compute SEs
  SE <- "Not calculated."
  covar <- "Not calculated."
  if (all(w[[1]] == w)){
    covar <- cppad_closed_estvar(Y, root, parts$offset, parts$Hess)
    SE <- sqrt(diag(covar))
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


