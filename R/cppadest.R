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
#' @param betaL NULL, a number, or a vector of length (p-1).
#' If NULL then the 1...(p-1) beta elements will be estimated.
#' If a number then the 1...(p-1) beta elements fixed at the given number.
#' If a vector, then the NA elements will be estimated and the others will be fixed at the supplied value.
#' @param betap NULL or a number.
#' If NULL then the pth element of beta will be estimated.
#' If a number, then the pth element of beta will be fixed at the given value.
#' @param control Control parameters to `Rcgmin`, of primary use is the tolerance parameter of the squared gradient size
#' @examples
#' model <- sec2_3model(1000)
#' estinfo <- ppi_cppad(model$sample, betap = -0.5, man = "Ralr", hsqfun = "ones")
#' misspecified <- ppi_cppad(model$sample, AL = "diag", bL = 0, betap = -0.5, man = "Ralr", hsqfun = "ones")
#' @export
ppi_cppad <- function(prop, AL = NULL, bL = NULL, Astar = NULL, betaL = NULL, betap = NULL,
                      pow = 1, man, hsqfun, acut = NULL, control = list(tol = 1E-20)){
  # process inputs
  stopifnot("matrix" %in% class(prop))
  p = ncol(prop)
  if ((man %in% c("Ralr", "Rclr", "Rmlr")) && (hsqfun != "ones")){
    warning("Manifold supplied has no boundary. Using hsqfun = 'ones' is strong recommended.")
  }
  stopifnot(pow == 1)

  theta <- ppi_cppad_thetaprocessor(p, AL, bL, Astar, betaL, betap)
  fixedtheta <- !is.na(theta)

  # prepare tapes
  pman <- pmanifold(man)
  if (man %in% c("Ralr", "Rclr", "Rmlr")){
    tapez <- rep(0.1, p - 1)
  } else {
    tapez <- rep(0.1, p)
  }
  thetatape <- theta  #must pass the fixed values as the taped value
  thetatape[!fixedtheta] <- 0.73 # any number will do!

  if (hsqfun == "ones"){acut = 1} #acut is igonred when hsqfun = ones
  pppi <- ptapell(tapez, thetatape,
                  llname = "ppi", pman,
                  fixedtheta = fixedtheta,
                  verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, p), rep(0.2, sum(!fixedtheta)),
                     pll = pppi, pman = pman,
                     hsqfun, acut = acut,
                     verbose = FALSE) #tape of the score function
  opt <- smest(smoppi, rep(0.2, sum(!fixedtheta)), prop, control = control)

  #process the theta and SE
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

ppi_cppad_thetaprocessor <- function(p, AL = NULL, bL = NULL, Astar = NULL, betaL = NULL, betap = NULL){
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
    } else if (is.na(AL)){
      ALprep = matrix(NA, nrow = p-1, ncol = p-1) #could also do nothing
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
    betapprep = betap
  }

  # combine above preparation into a vector, NA values to be estimated
  beta = c(betaLprep, betapprep)
  theta <- toPPIparamvec(ALprep, bLprep, beta)
  return(theta)
}
