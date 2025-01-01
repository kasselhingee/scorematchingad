#' @title Build CppAD Tapes for Score Matching
#' @family tape builders
#' @family generic score matching tools
#' @param thetatape_creator A function that generates tape values for theta. Must take a single argument, `n` the number for values to generate.
#' @param bdryw The name of the boundary weight function. "ones" for manifolds without boundary. For the simplex and positive orthant of the sphere, "prodsq" and "minsq" are possible - see [`ppi()`] for more information on these.
#' @param acut A parameter passed to the boundary weight function `bdryw`. Ignored for `bdryw = "ones"`.
#' @param verbose If `TRUE` more details are printed when taping. These details are for debugging and will likely be comprehensible only to users familiar with the source code of this package.
#' @description
#' For a parametric model family, the function `tape_smd()` generates `CppAD` tapes for the unnormalised log-density of the model family and of the score matching discrepancy function \eqn{A(z) + B(z) + C(z)} (defined in [`scorematchingtheory`]).
#' Three steps are performed by `tape_smd()`: first an object that specifies the manifold and any transformation to another manifold is created; then a tape of the unnormalised log-density is created; finally a tape of \eqn{A(z) + B(z) + C(z)} is created.
#' @details
#' To build a tape for the score matching discrepancy function, the `scorematchingad` first tapes the map from a point \eqn{z} on the `end` manifold to the value of the unnormalised log-density, where the independent variable is the \eqn{z}, the dynamic parameter is a vector of the parameters to estimate, and the remaining model parameters are fixed and not estimated.
#' This tape is then used to generate a tape for the score matching discrepancy function where the parameters to estimate are the independent variable.
#' 
#' The transforms of the manifold must be implemented in `C++` and selected by name.
#' @references \insertAllCited{}

#' @return
#' A list of:
#'   + an [`Rcpp_ADFun`] object containing a tape of the unnormalised log-density using the metric of the "`end`" manifold (that is the independent variable is on the `end` manifold).
#'   + an [`Rcpp_ADFun`] object containing a tape of the score matching discrepancy function with the non-fixed parameters of the model as the independent variable, and the measurements on the `end` manifold as the dynamic parameter.
#'   + some information about the tapes
#'
#' 
#' @examples
#' p <- 3
#' u <- rep(1/sqrt(p), p)
#' ltheta <- p #length of vMF parameter vector
#' intheta <- rep(NA, length.out = ltheta)
#' tapes <- tape_smd("sph", "identity", "sph", "vMF",
#'               ytape = u,
#'               usertheta = intheta,
#'               "ones", verbose = FALSE
#'               )
#' evaltape(tapes$lltape, u, runif(n = ltheta))
#' evaltape(tapes$smdtape, runif(n = ltheta), u)
#' 
#' u <- rep(1/3, 3)
#' tapes <- tape_smd("sim", "sqrt", "sph", "ppi",
#'               ytape = u,
#'               usertheta = ppi_paramvec(p = 3),
#'               bdryw = "minsq", acut = 0.01,
#'               verbose = FALSE
#'               )
#' evaltape(tapes$lltape, u, rppi_egmodel(1)$theta)
#' evaltape(tapes$smdtape, rppi_egmodel(1)$theta, u)
#' @export
tape_smd <- function(start, tran = "identity", end = start, ll,
                         ytape, usertheta,
                         bdryw = "ones", acut = 1,
                         thetatape_creator = function(n){seq(length.out = n)},
                         verbose = FALSE){

  if(all(!is.na(usertheta))){stop("All elements of theta are fixed")}

  tranman <- manifoldtransform(start, tran, end)

  # check bdryw and associated acut
  if (!(all(c(tran, end) == c("sqrt", "sph")) | all(c(tran, end) == c("identity", "sim")))){
    if (bdryw != "ones"){warning("Manifold supplied has no boundary. Using bdryw = 'ones' is strongly recommended.")}
  }
  if ((bdryw == "ones") && (abs(acut - 1) > 1E-8)){
    warning("The value of 'acut' is ignored for bdryw == 'ones'")
  }

  lltape <- tapell(ll = ll,
                    ytape = ytape,
                    usertheta = usertheta, 
                    thetatape_creator = thetatape_creator,
                    tranobj = tranman$tran)
  stopifnot(is.numeric(acut))
  smdtape <- tapesmd(uldtape = lltape,
                        tran = tranman$tran,
                        M = tranman$man,
                        weightname = bdryw,
                        acut = acut,
                        verbose = verbose)
  return(list(
    lltape = lltape,
    smdtape = smdtape,
    info = list(
      transform = tran,
      endmanifold = end,
      ulength = length(ytape),
      bdryw = bdryw,
      acut = acut
    )
  ))
}


