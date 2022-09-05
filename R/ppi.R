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
#' @param cW Specifies the tuning mutliplier `c` for computing Windham Weights. NULL for no robustness. Otherwise, easiest way specify `cW` is via [ppi_cW()] or [ppi_cW_auto()]. Use [ppi_paramvec()] greater customisation (at your own risk).
#' @examples
#' model <- ppi_egmodel(1000)
#' estinfo <- ppi(model$sample, paramvec = ppi_paramvec(betap = -0.5, p = ncol(model$sample)), trans = "alr", method = "cppad")
#' misspecified <- ppi(model$sample, paramvec = ppi_paramvec(bL = 0, betap = -0.5, p = ncol(model$sample)), trans = "alr", method = "direct")
#' @export
ppi <- function(Y, paramvec = NULL,
                pow = 1, trans, method = "direct", w = rep(1, nrow(Y)), cW = NULL,
                bdryweight = "ones", acut = NULL, #specific to some methods
                bdrythreshold = 1E-10, shiftsize = bdrythreshold, approxorder = 10, control = default_Rcgmin()#specific to cppad methods
                ){
  # process inputs
  stopifnot("matrix" %in% class(Y))
  p = ncol(Y)
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
  estimator <- function(...){stop("Estimator function needs to be defined")}
  fitfun <- NA

  controls <- splitcontrol(control)

  if (method == "direct"){
    if (man == "Ralr"){
      if (usertheta_ppi_alr_gengamma_compatible(usertheta)){
        firstfit <- ppi_alr_gengamma(Y, betap = tail(usertheta, 1), w = w) #any theta is fine
        fitfun <- "ppi_alr_gengamma"
        estimator <- function(Y, starttheta, isfixed, w){
           out <- ppi_alr_gengamma(Y, betap = tail(starttheta, 1), w = w)
           return(out$est$paramvec)
        }
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
          estimator <- function(Y, starttheta, isfixed, w){
             out <- as.vector(dir_sqrt_minimah(Y, acut = acut, w = w))
             return(t_sfi2u(out, starttheta, isfixed))
          }
        } else if (ppi_usertheta_estimator1_compatible_zerob(usertheta)){
          firstfit <- estimator1(Y,acut = acut,incb = 0,
                            beta = fromPPIparamvec(usertheta)$beta,
                            w= w, computeSE = TRUE)
          fitfun <- "estimator1_zerob"
          estimator <- function(Y, starttheta, isfixed, w){
             out <- estimator1(Y, acut = acut, incb = 0, beta = fromPPIparamvec(starttheta)$beta, w = w)
             return(out$est$paramvec)
          }
        } else if (ppi_usertheta_estimator1_compatible_incb(usertheta)){
          firstfit <- estimator1(Y,acut = acut,incb = 1,
                            beta = fromPPIparamvec(usertheta)$beta,
                            w= w, computeSE = TRUE)
          fitfun <- "estimator1_incb"
          estimator <- function(Y, starttheta, isfixed, w){
             out <- estimator1(Y, acut = acut, incb = 1, beta = fromPPIparamvec(starttheta)$beta, w = w)
             return(out$est$paramvec)
          }
        } else if (utheta_estimatorall1_betap_compatible(usertheta)){
          firstfit <- ppi_sqrt_minimah_full(Y, acut, tail(fromPPIparamvec(usertheta)$beta, 1),
                                            w)
          fitfun <- "ppi_sqrt_minimah_betap"
          estimator <- function(Y, starttheta, isfixed, w){
             out <- estimatorall1(Y, acut = acut,
                                  betap = tail(fromPPIparamvec(starttheta)$beta, 1),
                                  w = w)
             return(out$theta)
          }
        } else if (utheta_estimatorall1_full_compatible(usertheta)){
          firstfit <- ppi_sqrt_minimah_full(Y, acut, betap = NULL,
                                            w)
          fitfun <- "ppi_sqrt_minimah_full"
          estimator <- function(Y, starttheta, isfixed, w){
             out <- estimatorall1(Y, acut = acut,
                                  betap = NULL,
                                  w = w)
             return(out$theta)
          }
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
          estimator <- function(Y, starttheta, isfixed, w){
             out <- as.vector(dir_sqrt_prodh(Y, acut = acut, w = w))
             return(t_sfi2u(out, starttheta, isfixed))
          }
        } else if (ppi_usertheta_estimator1_compatible_zerob(usertheta)){
          firstfit <- ppi_sqrt_prodh_zerob(Y, acut, beta = fromPPIparamvec(usertheta)$beta, w)
          fitfun <- "ppi_sqrt_prodh_zerob"
          estimator <- function(Y, starttheta, isfixed, w){
             out <- estimator2(Y, acut = acut, incb = 0,
                               beta0 = fromPPIparamvec(starttheta)$beta,
                               w = w)
             return(out$theta)
          }
        } else if (ppi_usertheta_estimator1_compatible_incb(usertheta)){
          firstfit <- ppi_sqrt_prodh(Y, acut, beta = fromPPIparamvec(usertheta)$beta, w)
          fitfun <- "ppi_sqrt_prodh"
          estimator <- function(Y, starttheta, isfixed, w){
             out <- estimator2(Y, acut = acut, incb = 1,
                               beta0 = fromPPIparamvec(starttheta)$beta,
                               w = w)
             return(out$theta)
          }
        }
      }
    }
    if (is.na(fitfun)){
      warning("No direct estimator exists for parameter set. Using cppad.")
      method <- "cppad"
    }
  }
  if (method == "cppad"){
    stheta <- t_u2s_const(usertheta, 0.2)
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
    fitfun <- "cppad"
    estimator <- function(Y, starttheta, isfixed, w){
      out <- ppi_cppad(Y, stheta = stheta, isfixed = isfixed,
               bdrythreshold = bdrythreshold,
               shiftsize = shiftsize,
               approxorder = approxorder,
               pow = pow,
               man = man,
               weightname = bdryweight,
               acut = acut,
               control = controls$Rcgmin,
               w = w)
      return(out$theta)
    }
  }

  #### No Robustness, return first fit ####
  if (is.null(cW)){
     firstfit$info$method <- fitfun
     return(firstfit)
  }

  #### Do Windham Robustness ####
  stopifnot(length(cW) == length(usertheta))
  stopifnot(is.numeric(cW))
  stopifnot(all((cW * usertheta)[t_u2i(usertheta)] == 0)) #all elements of cW for fixed non-zero parameters should be zero.

  ldenfun <- function(Y, theta){ #here theta is the usual parameters of PPI model from
    mats <- fromPPIparamvec(theta, p = ncol(Y))
    return(drop(dppi(Y, beta0=mats$beta, ALs = mats$ALs, bL = mats$bL)))
  }

  est <- windham_raw(prop = Y,
                     cW = cW,
                     ldenfun = ldenfun,
                     estimatorfun = estimator,
                     starttheta = firstfit$est$paramvec,
                     isfixed = t_u2i(usertheta),
                     originalcorrectionmethod = FALSE, #for variable cW
                     fpcontrol = controls$fp)

  return(c(est, fitfun = fitfun))
}

