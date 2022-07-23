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
#' @param hsqfun The h-squared function for down weighting as measurements approach the manifold boundary.
#' @param acut The threshold in `hsqfun` to avoid over-weighting measurements interior to the simplex
#' @param AL NULL, a p-1 x p-1 symmetric matrix, a number, or "diag".
#' If NULL then AL will be estimated.
#' If a matrix, then the NA elements will be estimated and the others will be fixed at the supplied value (i.e. not estimated).
#' If a single number, then AL will be fixed as a matrix of the given value.
#' If "diag" then the non-diagonal elements of AL will be fixed to 0.
#' @param bL NULL, a number, or a vector of length (p-1).
#' If NULL, then bL will be estimated.
#' If a number, then bL will be fixed at the supplied value.
#' If a vector, then the NA elements will be estimated and the others will be fixed at the supplied value.
#' @param Astar  NULL or a p x p matrix.
#' If non-null, then overrides AL and bL.
#' If a matrix, all elements must be provided and Astar will be fixed in the estimation
#' (This is because transforming to AL and bL from an incomplete Astar appears impossible).
#' @param beta NULL, a number, or a vector of length p. If non-null then overrides `betaL` and `betap` arguments.
#' If NULL then the elements of the beta vector will be estimated.
#' If a number then the beta elements will be fixed at the given number.
#' If a vector, then the NA elements will be estimated and the others will be fixed at the supplied value.
#' @param betaL NULL, a number, or a vector of length (p-1). If non-null then overrides `beta` argument.
#' If NULL then the 1...(p-1) beta elements will be estimated.
#' If a number then the 1...(p-1) beta elements fixed at the given number.
#' If a vector, then the NA elements will be estimated and the others will be fixed at the supplied value.
#' @param betap NULL or a number. If non-null then overrides `beta` argument.
#' If NULL then the pth element of beta will be estimated.
#' If a number, then the pth element of beta will be fixed at the given value.
#' @param control Control parameters to `Rcgmin`, of primary use is the tolerance parameter of the squared gradient size
#' @param bdrythreshold For measurements close to the boundary of the simplex Taylor approximation is applied.
#' @param shiftsize Measurements close to the boundary are shifted by this distance for Taylor approximation.
#' @param approxorder Order of the Taylor approximation
#' @param method `direct` for estimates calculated directly where possible (*list them*) or `cppad` to find the score matching estimates using automatic differentiation and the `Rcgmin()` iterative solver.
#' @examples
#' model <- ppi_egmodel(1000)
#' estinfo <- ppi(model$sample, betap = -0.5, man = "Ralr", weightname = "ones")
#' misspecified <- ppi(model$sample, AL = "diag", bL = 0, betap = -0.5, man = "Ralr", weightname = "ones")
#' @export
ppi <- function(Y, AL = NULL, bL = NULL, Astar = NULL, beta = NULL, betaL = NULL, betap = NULL,
                pow = 1, trans, method = "direct", w = rep(1, nrow(Y)),
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

  usertheta <- ppi_cppad_thetaprocessor(p, AL, bL, Astar, beta, betaL, betap)
  out <- list()
  fitfun <- NA

  if (method == "direct"){
    if (man == "Ralr"){
        if (usertheta_estimatorlog_weight_compatible(usertheta)){
        out <- estimatorlog_weight(Y, betap = usertheta[length(usertheta)], weightW = w) #any theta is fine
        fitfun <- "estimatorlog_weight"
        }
    }
    if (man == "sphere"){ # a number of methods implemented
      if (bdryweight == "minsq"){
        if (ppi_usertheta_for_dir_sqrt_minimah(usertheta)){
          out$est <- dir_sqrt_minimah(Y, acut = acut, w = w)
          fitfun <- "dir_sqrt_minimah"
        } else if (ppi_usertheta_estimator1_compatible_zerob(usertheta)){
          out <- estimator1(Y,acut = acut,incb = 0,
                            beta0 = fromPPIparamvec(usertheta)$beta,
                            w= w)
          fitfun <- "estimator1_zerob"
        } else if (ppi_usertheta_estimator1_compatible_incb(usertheta)){
          out <- estimator1(Y,acut = acut,incb = 1,
                            beta0 = fromPPIparamvec(usertheta)$beta,
                            w= w)
          fitfun <- "estimator1_incb"
        } else if (utheta_estimatorall1_betap_compatible(usertheta)){
          out <- estimatorall1(Y, acut = acut,
                            betap = tail(fromPPIparamvec(usertheta)$beta, 1),
                            w= w)
          fitfun <- "estimatorall1_betap"
        } else if (utheta_estimatorall1_full_compatible(usertheta)){
          out <- estimatorall1(Y, acut = acut,
                               betap = NULL,
                               w= w)
          fitfun <- "estimatorall1_full"
        }
      }

      if (bdryweight == "prodsq"){
        if (ppi_usertheta_for_dir_sqrt_minimah(usertheta)){
          out$est <- dir_sqrt_prodh(Y, acut = acut, w = w)
          fitfun <- "dir_sqrt_prodh"
        } else if (ppi_usertheta_estimator1_compatible_zerob(usertheta)){
          out <- estimator2(Y,acut = acut,incb = 0,
                            beta0 = fromPPIparamvec(usertheta)$beta,
                            w= w)
          fitfun <- "estimator2_zerob"
        } else if (ppi_usertheta_estimator1_compatible_incb(usertheta)){
          out <- estimator2(Y,acut = acut,incb = 1,
                            beta0 = fromPPIparamvec(usertheta)$beta,
                            w= w)
          fitfun <- "estimator2_incb"
        }
      }
    }
    if (is.na(fitfun)){
      warning("No direct estimator exists for parameter set. Using cppad.")
      method <- "cppad"
    }
  }
  if (method == "cppad"){
    out <- ppi_cppad(Y, usertheta, bdrythreshold, shiftsize, approxorder, pow, man, bdryweight, acut, control)
    fitfun <- "cppad"
  }

  return(c(fitfun = fitfun,
           out))
}


ppi_cppad <- function(prop, usertheta,
                      bdrythreshold = 1E-10, shiftsize = bdrythreshold, approxorder = 10,
                      pow = 1, man, weightname = hsqfun, acut = NULL, control = default_Rcgmin(), hsqfun = NULL){
  # process inputs
  stopifnot("matrix" %in% class(prop))
  p = ncol(prop)

  theta <- usertheta

  if (!(man %in% c("simplex", "sphere"))){
    if (weightname != "ones"){warning("Manifold supplied has no boundary. Setting weightname to 'ones'.")}
  }
  if (weightname == "ones"){
    if (!is.null(acut)){warning("The value of 'acut' is ignored for weightname == 'ones'")}
    acut <- 1 #set just for passing to Cpp
  }

  # prepare tapes
  tapes <- buildsmotape(man, "ppi",
                rep(1/p, p), theta,
                weightname = weightname,
                acut = acut, verbose = FALSE
                )

  # split data into boundary and interior
  datasplit <- simplex_boundarysplit(prop, bdrythreshold = bdrythreshold, shiftsize = shiftsize)

  opt <- smest(tapes$smotape, rep(0.2, sum(is.na(theta))), datasplit$interior,
               uboundary = datasplit$uboundary, boundaryapprox = datasplit$boundaryapprox,
               approxorder = approxorder,
               control = control)

  #process the theta and SE
  fixedtheta <- !is.na(theta)
  thetaest <- theta
  thetaest[!fixedtheta] <- opt$par
  SE <- fixedtheta * 0
  SE[!fixedtheta] <- opt$SE

  # make output
  list(
    prop = prop,
    est = c(list(theta = thetaest),
            fromPPIparamvec(thetaest, p)),
    SE = c(list(theta = SE),
           fromPPIparamvec(SE, p)),
    smval = opt$value,
    sqgradsize = opt$sqgradsize,
    counts = opt$counts,
    convergence = opt$convergence,
    message = opt$convergence
    )
}

ppi_cppad_thetaprocessor <- function(p, AL = NULL, bL = NULL, Astar = NULL, beta = NULL, betaL = NULL, betap = NULL){
  # initialise parameter objects
  bLprep = rep(NA, p-1)
  betaLprep = rep(NA, p-1)
  betapprep = NA
  # AL, bL and A
  if (!is.null(Astar)){
    if (!(is.null(AL) & is.null(bL))){warning("AL, bL and Astar supplied. Astar argument will override AL and bL.")}
    translated <- fromAstar(Astar)
    ALprep <- translated$AL
    bLprep <- translated$bL
  } else {
    # AL
    if (is.null(AL)){
      ALprep = matrix(NA, nrow = p-1, ncol = p-1) #could also do nothing
    } else if (is.matrix(AL)){
    # If a matrix, then the NA elements will be estimated and the others will be fixed at the supplied value (i.e. not estimated).
      if(!isSymmetric.matrix(AL)){stop("AL must be symmetric.")}
      ALprep = AL
    } else if (is.numeric(AL)){#' If a single number, then AL will be fixed as a matrix of the given value.
      ALprep = matrix(AL, nrow = p-1, ncol = p-1)
    } else if (is.character(AL)){#' If "diag" then the non-diagonal elements of AL will be fixed to 0.
      stopifnot((AL == "diag") | (AL == "d") | (AL == "diagonal"))
      ALprep <- matrix(0, nrow = p-1, ncol = p-1)
      diag(ALprep) <- NA
    } else if (is.logical(AL)){
      ALprep = matrix(AL, nrow = p-1, ncol = p-1) #covers NA, TRUE, and FALSE
    } else {
      stop("AL is not of required type.")
    }
    #bL
    # If a number, then bL will be fixed at the supplied value.
    if (!is.null(bL)){
      if (!is.vector(bL, mode = "any")){stop("bL must be a vector or value")}
      if (length(bL) == 1){
        bLprep = rep(bL, p-1)
      } else {#' If a vector, then the NA elements will be estimated and the others will be fixed at the supplied value.
        if (length(bL) != p-1){stop("bL must have length p-1")}
        bLprep = bL
      }
    }
  }

  # beta
  if (!is.null(betaL)){
    stopifnot(is.vector(betaL, "numeric") | is.vector(betaL, "logical"))
    stopifnot(is.null(beta))
    if (length(betaL) == 1){
      # If a number then the 1...(p-1) beta elements fixed at the given number.
      betaLprep = rep(betaL, p-1)
    } else if (length(betaL) == p-1){
      # If a vector, then the NA elements will be estimated and the others will be fixed at the supplied value.
      betaLprep = betaL
    } else {
      stop("betaL must have length p-1")
    }
  }
  if (!is.null(betap)){
    stopifnot(length(betap) == 1)
    stopifnot(is.numeric(betap) | is.logical(betap))
    stopifnot(is.null(beta))
    betapprep = betap
  }
  if (!is.null(beta)){
    stopifnot(is.null(betaL))
    stopifnot(is.null(betap))
    if (is.matrix(beta)){beta <- drop(beta)}
    stopifnot(is.vector(beta, "numeric") | is.vector(beta, "logical"))
    if (length(beta) == 1){
      betaLprep = rep(beta, p-1)
      betapprep = beta
    } else if (length(beta) == p){
      # If a vector, then the NA elements will be estimated and the others will be fixed at the supplied value.
      betaLprep = beta[1:(p-1)]
      betapprep = beta[p]
    } else {
      stop("beta must have length p")
    }
  }

  # combine above preparation into a vector, NA values to be estimated
  beta = c(betaLprep, betapprep)
  theta <- toPPIparamvec(ALprep, bLprep, beta)
  return(theta)
}
