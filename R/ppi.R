#' @title Score-Matching Estimation of PPI Parameters
#' @description For certain situations computes the score matching estimate directly (e.g. \insertCite{scealy2022sc}{cdabyppi}), otherwise iteratively minimises the *Hyvarinen divergence* \insertCite{@Equation 2, @hyvarinen2005es}{cdabyppi} using derivatives computed by CppAD and [Rcgmin::Rcgmin()].

#' @details
#' Estimation may be performed via transformation onto Euclidean space, the positive quadrant of the sphere, or without any transformation. In the latter two situations there is a boundary and *weighted Hyvarinen divergence* \insertCite{@Equation 7, @scealy2022sc}{cdabyppi} is used.
#'
#' Direct estimates are available for the following situations
#' + `trans='alr'` and `betap` supplied (and typically positive)
#' + `trans='sqrt'` and ....

#' @param trans The name of the transformation: 'alr' (additive log ratio), 'sqrt' or 'none'.
#' @param pow The power of `u` in the PPI density - by default `pow` is `1`. NOT YET IMPLEMENTED
#' @param Y A matrix of measurements. Each row is a measurement, each component is a dimension of the measurement.
#' @param paramvec Optionally a standard-form vector of the PPI models parameters. Non-NA values are fixed, NA-valued elements are estimated. Generate `paramvec` easily using [ppi_paramvec()].
#' @param hsqfun The h-squared function for down weighting as measurements approach the manifold boundary.
#' @param acut The threshold in `hsqfun` to avoid over-weighting measurements interior to the simplex
#' @param control Control parameters to `Rcgmin`, of primary use is the tolerance parameter of the squared gradient size
#' @param bdrythreshold For measurements close to the boundary of the simplex Taylor approximation is applied.
#' @param shiftsize Measurements close to the boundary are shifted by this distance for Taylor approximation.
#' @param approxorder Order of the Taylor approximation
#' @param method `direct` for estimates calculated directly where possible (*list them*) or `cppad` to find the score matching estimates using automatic differentiation and the `Rcgmin()` iterative solver.
#' @param paramvec_start Only for method `cppad`. The starting guess for the iterative solver with possibly NA values for the fixed (not-estimated) elements. Generate `paramvec` easily using [ppi_paramvec()].
#' @examples
#' model <- ppi_egmodel(1000)
#' estinfo <- ppi(model$sample, paramvec = ppi_paramvec(betap = -0.5, p = ncol(model$sample)), trans = "alr", method = "cppad")
#' misspecified <- ppi(model$sample, paramvec = ppi_paramvec(bL = 0, betap = -0.5, p = ncol(model$sample)), trans = "alr", method = "direct")
#' @export
ppi <- function(Y, paramvec = NULL,
                pow = 1, trans, method = "direct", w = rep(1, nrow(Y)),
                bdryweight = "ones", acut = NULL, #specific to some methods
                bdrythreshold = 1E-10, shiftsize = bdrythreshold, approxorder = 10, control = default_Rcgmin(), paramvec_start = NULL#specific to cppad methods
                ){
  # process inputs
  stopifnot("matrix" %in% class(Y))
  p = ncol(Y)
  if (any(Y < 0)){
     warning(sprintf("Y contains %i negative values.", sum(Y < 0)))
  }
  if (any(rowSums(Y) != 1)){
     warning(paste("Y contains measurement that don't add to 1. Largest discrepancy is", max(abs(rowSums(Y) - 1))))
  }

  stopifnot(pow == 1)
  stopifnot(trans %in% c("alr", "sqrt", "none"))
  man <- switch(trans,
           alr = "Ralr",
           sqrt = "sphere",
           none = "simplex")

  if (is.null(paramvec)){usertheta <- rep(NA, ppithetalength(p))}
  else {usertheta <- paramvec}

  stopifnot(length(usertheta) == ppithetalength(p))

  firstfit <- list()
  fitfun <- NA

  controls <- splitcontrol(control)

  # Switch between the different methods
  if (method == "direct"){
    if (man == "Ralr"){
      if (usertheta_ppi_alr_gengamma_compatible(usertheta)){
        firstfit <- ppi_alr_gengamma(Y, betap = tail(usertheta, 1), w = w) #any theta is fine
        fitfun <- "ppi_alr_gengamma"
      }
    }
    if (man == "sphere"){ # a number of methods implemented
      if (bdryweight == "minsq"){
        if (ppi_usertheta_for_dir_sqrt_minimah(usertheta)){
          betaest <- as.vector(dir_sqrt_minimah(Y, acut = acut, w = w))
          estparamvec <- t_fu2t(betaest, usertheta)
          firstfit$est <- c(list(paramvec = estparamvec),
                            fromPPIparamvec(estparamvec))
          firstfit$SE <- "Not calculated."
          fitfun <- "dir_sqrt_minimah"
        } else if (ppi_usertheta_estimator1_compatible_zerob(usertheta)){
          firstfit <- estimator1(Y,acut = acut,incb = 0,
                            beta = fromPPIparamvec(usertheta)$beta,
                            w= w, computeSE = TRUE)
          fitfun <- "estimator1_zerob"
        } else if (ppi_usertheta_estimator1_compatible_incb(usertheta)){
          firstfit <- estimator1(Y,acut = acut,incb = 1,
                            beta = fromPPIparamvec(usertheta)$beta,
                            w= w, computeSE = TRUE)
          fitfun <- "estimator1_incb"
        } else if (utheta_estimatorall1_betap_compatible(usertheta)){
          firstfit <- ppi_sqrt_minimah_full(Y, acut, tail(fromPPIparamvec(usertheta)$beta, 1),
                                            w)
          fitfun <- "ppi_sqrt_minimah_betap"
        } else if (utheta_estimatorall1_full_compatible(usertheta)){
          firstfit <- ppi_sqrt_minimah_full(Y, acut, betap = NULL,
                                            w)
          fitfun <- "ppi_sqrt_minimah_full"
        }
      }

      if (bdryweight == "prodsq"){
        if (ppi_usertheta_for_dir_sqrt_minimah(usertheta)){
          betaest <- dir_sqrt_prodh(Y, acut = acut, w = w)
          fitfun <- "dir_sqrt_prodh"
          estparamvec <- t_fu2t(betaest, usertheta)
          firstfit$est <- c(list(paramvec = estparamvec),
                            fromPPIparamvec(estparamvec))
          firstfit$SE <- "Not calculated."
        } else if (ppi_usertheta_estimator1_compatible_zerob(usertheta)){
          firstfit <- ppi_sqrt_prodh_zerob(Y, acut, beta = fromPPIparamvec(usertheta)$beta, w)
          fitfun <- "ppi_sqrt_prodh_zerob"
        } else if (ppi_usertheta_estimator1_compatible_incb(usertheta)){
          firstfit <- ppi_sqrt_prodh(Y, acut, beta = fromPPIparamvec(usertheta)$beta, w)
          fitfun <- "ppi_sqrt_prodh"
        }
      }
    }
    if (is.na(fitfun)){
      warning("No direct estimator exists for parameter set. Using cppad.")
      method <- "cppad"
    }
  }
  if (method == "cppad"){
    if (is.null(paramvec_start)){stheta <- t_u2s_const(usertheta, 0.2)}
    else {stheta <- t_us2s(usertheta, paramvec_start)}
    isfixed <- t_u2i(usertheta)
    firstfit <- ppi_cppad(Y, stheta = stheta, isfixed = isfixed,
               bdrythreshold = bdrythreshold,
               shiftsize = shiftsize,
               approxorder = approxorder,
               pow = pow,
               man = man,
               weightname = bdryweight,
               acut = acut,
               control = controls$Rcgmin,
               w = w)
    #refactor results to fit with ppi() standard output
    firstfit$est <- c(list(paramvec = firstfit$est$theta),
                      fromPPIparamvec(firstfit$est$theta))
    firstfit$SE <- c(list(paramvec = firstfit$SE$theta),
                      fromPPIparamvec(firstfit$SE$theta))
    names(firstfit)
    firstfit$info <- firstfit[setdiff(names(firstfit), c("theta", "prop", "est", "SE", "info"))]
    firstfit[setdiff(names(firstfit), c("est", "SE", "info"))] <- NULL
    #
    fitfun <- "cppad"
  }

  firstfit$info$method <- fitfun
  return(firstfit)
}

