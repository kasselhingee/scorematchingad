#' @title Estimation of PPI Parameters via CppAD
#' @description Uses CppAD to compute the score matching objective value and `Rcgmin` to minimise it.
#' Measurements can be automatically transformed to a number of different manifolds (using different transformations).
#' Different weight functions are also available.
#' @param man Manifold name (includes tranformation)
#' @param pow The power of `u` in the PPI density - by default `pow` is `1`. NOT YET IMPLEMENTED
#' @param prop A matrix of measurements. Each row is a measurement, each component is a dimension of the measurement.
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
#' model <- sec2_3model(1000)
#' estinfo <- ppi_cppad(model$sample, betap = -0.5, man = "Ralr", weightname = "ones")
#' misspecified <- ppi_cppad(model$sample, AL = "diag", bL = 0, betap = -0.5, man = "Ralr", weightname = "ones")
#' @export
ppi <- function(Y, AL = NULL, bL = NULL, Astar = NULL, beta = NULL, betaL = NULL, betap = NULL,
                pow = 1, man, method = "direct", w = NULL,
                bdryweight = "ones", acut = NULL, #specific to some methods
                bdrythreshold = 1E-10, shiftsize = bdrythreshold, approxorder = 10, control = default_Rcgmin()#specific to cppad methods
                ){
  # process inputs
  stopifnot("matrix" %in% class(Y))
  p = ncol(Y)
  stopifnot(pow == 1)

  usertheta <- ppi_cppad_thetaprocessor(p, AL, bL, Astar, beta, betaL, betap)
  if (method == "cppad"){
    out <- ppi_cppad(Y, usertheta, bdrythreshold, shiftsize, approxorder, pow, man, bdryweight, acut, control)
  }

  return(out)
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
    #' If a matrix, then the NA elements will be estimated and the others will be fixed at the supplied value (i.e. not estimated).
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
    #' If a number, then bL will be fixed at the supplied value.
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
      #' If a number then the 1...(p-1) beta elements fixed at the given number.
      betaLprep = rep(betaL, p-1)
    } else if (length(betaL) == p-1){
      #' If a vector, then the NA elements will be estimated and the others will be fixed at the supplied value.
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
    stopifnot(is.vector(beta, "numeric") | is.vector(beta, "logical"))
    if (length(beta) == 1){
      betaLprep = rep(beta, p-1)
      betapprep = beta
    } else if (length(beta) == p){
      #' If a vector, then the NA elements will be estimated and the others will be fixed at the supplied value.
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
