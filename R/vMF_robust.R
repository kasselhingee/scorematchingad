#' @title Robust fitting of von Mises Fisher
#' @description
#' Robust fitting for von Mises Fisher distribution using [`Windham()`], 
#' for the full parameter vector \eqn{\kappa\mu} with `vMF_robust()` or just the concentration with `vMF_kappa_robust()` assuming a mean direction of `c(1, 0, 0, ....)`.
#' @param cW Tuning constants for each parameter in the vMF parameter vector. If a single number then the constant is the same for each element of the parameter vector.
#' @param Y A matrix of observations in Cartesian coordinates.
#' @param ... Passed to [`Windham()`] and then passed onto either [`vMF()`] or [`vMF_kappa()`].
#' @export
vMF_robust <- function(Y, cW, ...){
  ellipsis::check_dots_used()
  extraargs <- list(...)

  # user friendly cW
  if (length(cW) == 1){
    if (!is.null(extraargs$paramvec)){
      isfixed = t_u2i(extraargs$paramvec)
      cW <- cW * !isfixed
    } else {
      cW <- rep(cW, ncol(Y))
    }
  }


  ldenfun <- function(Y, theta){ #here theta is km
    k <- sqrt(sum(theta^2))
    m <- theta/k
    return(drop(Directional::dvmf(Y, k, m, logden = TRUE)))
  }
  est <- WindhamRobust(Y = Y,
                     estimator = vMF,
                     ldenfun = ldenfun,
                     cW = cW,
                     ...)
  out <- list(
    est = c(list(paramvec = est$theta), vMF_fromparamvec(est$theta)),
    SE = "Not calculated.",
    info = est$optim
  )
  return(out)
}

#' @rdname 
#' @export
vMF_kappa_robust <- function(Y, cW, ...){
  extraargs <- list(...)
  Y <- vMF_stdY(Y, w = extraargs$w)
  ellipsis::check_dots_used()
  ldenfun <- function(Y, theta){ #here theta is k and m is c(1, 0, ...)
    k <- theta
    m <- c(1 , rep(0, ncol(Y) - 1))
    return(drop(Directional::dvmf(Y, k, m, logden = TRUE)))
  }
  est <- WindhamRobust(Y = Y,
                     estimator = vMF_kappa,
                     ldenfun = ldenfun,
                     cW = cW,
                     ...)
  out <- list(
    est = list(paramvec = est$theta, k = est$theta),
    SE = "Not calculated.",
    info = est$optim
  )
  return(out)
}
