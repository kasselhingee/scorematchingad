#' @title Compute score matching discrepancy value, gradient, and Hessian for a PPI Model
#' @family PPI model tools
#' @description Using similar arguments to [`ppi()`], compute values related to score matching. See [`smvalues_tape()`]. The gradient offset is also computed (see [`quadratictape_parts()`]. 
#' @inheritParams ppi
#' @param evalparam The parameter set to evaluate the score matching values.
#' This is different to `paramvec`, which *fixes* parameters and constrains the estimation.
#'  All elements of `evalparam` must be non-NA, and if there are any parameter fixed by `paramvec` then `evalparam` must match them (and it will warn if not). 
#' @param average If TRUE return the (weighted average) of the measurements, otherwise return the values for each measurement.
#' @return
#' A list of 
#'  + `obj` the score matching discrepancy value
#'  + `grad` the gradient of the score matching discrepancy
#'  + `hess` the Hessian of the score matching discrepancy
#'  + `offset` gradient offset (see [`quadratictape_parts()`])
#' @export
ppi_smvalues <- function(Y, paramvec = NULL, evalparam,
                trans, method = "closed", w = rep(1, nrow(Y)),
                divweight = "ones", acut = NULL, #specific to some methods
                bdrythreshold = 1E-10, shiftsize = bdrythreshold, approxorder = 10, #specific to cppad methods
                average = TRUE
                ){
  ###### process inputs #####
  stopifnot("matrix" %in% class(Y))
  p = ncol(Y)
  numneg <- sum(Y < 0)
  if (numneg > 0){
     warning(sprintf("Y contains %i negative values.", numneg))
  }
  sum_m1 <- max(abs(rowSums(Y) - 1))
  if (sum_m1 > 1E-15){
     warning(sprintf("Y contains measurement that don't add to 1. Largest discrepancy is %s.", sum_m1))
  }

  stopifnot(trans %in% c("alr", "sqrt", "clr", "none"))
  man <- switch(trans,
           alr = "Euc",
           clr = "Hn111",
           sqrt = "sph",
           none = "sim")
  if (trans == "none"){trans <- "identity"}

  if (man %in% c("sim", "sph")){
    if (divweight == "ones"){stop("Manifold supplied has a boundary - set divweight to something that isn't 'ones'")}
  } else {
    if (divweight != "ones"){warning("Manifold supplied has no boundary. Setting divweight to 'ones'.")}
  }
  if (divweight == "ones"){
    if (!is.null(acut)){warning("The value of 'acut' is ignored for divweight == 'ones'")}
    acut <- 1 #set just for passing to CppAD
  }

  if (is.null(paramvec)){paramvec <- rep(NA, ppithetalength(ncol(Y)))}

  tapes <- buildsmotape(
     start = "sim",
     tran = trans,
     man = man,
     llname = "ppi",
     ytape =  rep(1/p, p),
     usertheta = paramvec,
     divweight = divweight,
     acut = acut,
     verbose = FALSE)
  smotape <- tapes$smotape

  # find boundary points and their approximation centres
  isbdry <- simplex_isboundary(Y, bdrythreshold)
  Yapproxcentres <- Y 
  Yapproxcentres[!isbdry, ] <- NA
  Yapproxcentres[isbdry, ] <- simplex_boundaryshift(Y[isbdry, , drop = FALSE], shiftsize = shiftsize)

  # gradient values
  if (average){
    valgradhess <- smvalues_tape_wsum(smotape, xmat = Y, pmat = t_ut2f(paramvec, evalparam), 
                              xcentres = Yapproxcentres,
                              w = w,
                              approxorder = approxorder)
  } else {
    valgradhess <- smvalues_tape(smotape, xmat = Y, pmat = t_ut2f(paramvec, evalparam), 
                                      xcentres = Yapproxcentres,
                                      approxorder = approxorder)
  }

  # quadratic simplification of divergence discrepancy
  quadparts <- quadratictape_parts(smotape, tmat = Y, tcentres = Yapproxcentres, approxorder = approxorder)

  # combine:
  if (average){
    quadparts_wsum <- lapply(quadparts, wcolSums, w = w)
  
    if (is.null(w)){
      normaliser <- nrow(Y)
    } else {
      normaliser <- sum(w)
    }
    out <- lapply(c(valgradhess, quadparts_wsum["offset"]), function(x){x/normaliser})
  } else {
    out <- c(valgradhess, quadparts["offset"])
  }
  return(out)
}
