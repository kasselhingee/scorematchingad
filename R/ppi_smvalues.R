#' @title Compute score matching objective value, gradient, and Hessian for a PPI Model
#' @family PPI model tools
#' @description Using similar arguments to [`ppi()`], compute values related to score matching. See [`tape_smvalues()`]. The gradient offset is also computed (see [`quadratictape_parts()`]. 
#' @inheritParams ppi
#' @param evalparam The parameter set to evaluate the score matching values.
#' This is different to `paramvec`, which *fixes* parameters and constrains the estimation.
#'  All elements of `evalparam` must be non-NA, and if there are any parameter fixed by `paramvec` then `evalparam` must match them (and it will warn if not). 
#' @return
#' A list of 
#'  + `obj` the score matching objective value
#'  + `grad` the gradient of the score matching objective
#'  + `hess` the Hessian of the score matching objective
#'  + `offset` gradient offset (see [`quadratictape_parts()`])
#' @export
ppi_smvalues <- function(Y, paramvec = NULL, evalparam,
                trans, method = "closed", w = rep(1, nrow(Y)),
                divweight = "ones", acut = NULL, #specific to some methods
                bdrythreshold = 1E-10, shiftsize = bdrythreshold, approxorder = 10 #specific to cppad methods
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
           alr = "Ralr",
           clr = "Hclr",
           sqrt = "sphere",
           none = "simplex")
  if (!(man %in% c("simplex", "sphere"))){
    if (divweight != "ones"){warning("Manifold supplied has no boundary. Setting divweight to 'ones'.")}
  }
  if (divweight == "ones"){
    if (!is.null(acut)){warning("The value of 'acut' is ignored for divweight == 'ones'")}
    acut <- 1 #set just for passing to CppAD
  }

  if (is.null(paramvec)){paramvec <- rep(NA, ppithetalength(ncol(Y)))}

  tapes <- buildsmotape(
     manifoldname = man,
     llname = "ppi",
     ytape =  rep(1/p, p),
     usertheta = paramvec,
     weightname = divweight,
     acut = acut,
     verbose = FALSE)
  smotape <- tapes$smotape

  # find boundary points and their approximation centres
  isbdry <- simplex_isboundary(Y, bdrythreshold)
  Yapproxcentres <- Y 
  Yapproxcentres[!isbdry, ] <- NA
  Yapproxcentres[isbdry, ] <- simplex_boundaryshift(Y[isbdry, , drop = FALSE], shiftsize = shiftsize)

  valgradhess <- tape_smvalues_wsum(smotape, xmat = Y, pmat = t_ut2f(paramvec, evalparam), 
                            xcentres = Yapproxcentres,
                            w = w,
                            approxorder = approxorder)
  quadparts <- quadratictape_parts(smotape, tmat = Y, tcentres = Yapproxcentres, approxorder = approxorder)
  quadparts_wsum <- lapply(quadparts, wcolSums, w = w)

  if (is.null(w)){
    normaliser <- nrow(Y)
  } else {
  normaliser <- sum(w)
  }
  out <- lapply(c(valgradhess, quadparts_wsum["offset"]), function(x){x/normaliser})
  return(out)
}
