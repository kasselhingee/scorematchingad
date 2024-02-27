#' @param tran A transform object (of type `Rcpp_transform_ad`), typically created by [`manifoldtransform()`].
#' @param llname The name of the log-likelihood function to tape
#' @param ytape An example measurement value to use for creating the tape. In the natural manifold of the log-likelihood function. `ytape` will be converted to the manifold according to the `toM()` method for `tran` before taping. 
#' Please ensure that `ytape` is the interior of the manifold, and it is probably best if all components of `tran$toM(ytape)` are non-zero.
#' @param usertheta A vector of parameter elements for the likelihood function. `NA` elements will become *dynamic parameters*. Other elements will be fixed at the provided value. The length of `usertheta` must be the correct length for the log-likelihood, no checking is conducted.
#' @param thetatape_creator A function that accepts an integer `n`, and returns a vector of `n` length. The function is used to fill in the `NA` elements of `usertheta` when building the tapes. Please ensure that the values filled by `thetatape_creator` lead to plausible parameter vectors for the chosen log-likelihood.
#' @describeIn buildsmotape Creates a `CppAD` tape of the log-likelihood as a function of values on the provided manifold. The log-likelihood accounts for the change of measure between the natural manifold of the log-likelihood and the manifold with transformation provided by `tran`.
#' @return 
#' `tapell()` returns an [`ADFun`] object with two additional attributes accessed via `attr()`:  
#'  + `ytape` The value of `ytape`
#'  + `tran` The name of the transform specified in `tran`.
#' @examples 
#' maninfo <- manifoldtransform("sim", "sqrt", "sph")
#' ppitape <- tapell(llname = "ppi",
#'                   ytape = c(0.2, 0.3, 0.5),
#'                   usertheta = ppi_paramvec(p = 3), 
#'                   tran = maninfo$tran) 
#' pForward0(ppitape$ptr, 
#'   sqrt(rep(1/3, 3)), 
#'   ppi_paramvec(p = 3, AL=0, bL=0, beta=c(0,0,0.5)))
#' pJacobian(ppitape$ptr, 
#'   sqrt(rep(1/3, 3)), 
#'   ppi_paramvec(p = 3, AL=0, bL=0, beta=c(0,0,0.5)))
#' pHessian(ppitape$ptr, 
#'   sqrt(rep(1/3, 3)), 
#'   ppi_paramvec(p = 3, AL=0, bL=0, beta=c(0,0,0.5)))
#' @section Warning: There is limited checking of the inputs.
#' @export
tapell <- function(llname,
                   ytape,
                   usertheta,
                   tran,
                   thetatape_creator = function(n){seq(length.out = n)},
                   verbose = FALSE){
  stopifnot(inherits(tran, "Rcpp_transform_ad"))
  starttheta <- t_u2s(usertheta, filler = thetatape_creator)
  ztape <- tran$toM(ytape) #the value of ytape transformed to the manifold
  lltape <- ptapell(ztape, starttheta,
                    llname = llname, 
                    tran = tran,
                    fixedtheta = t_u2i(usertheta),
                    verbose = verbose)
  out <- ADFun$new(ptr = lltape,
                   name = paste(tran$name(), llname, sep = "-"),
                   xtape = ztape,
                   dyntape =  as.numeric(starttheta[!t_u2i(usertheta)]),
                   usertheta = as.numeric(usertheta))

  attr(out, "ytape") <- ytape
  attr(out, "tran") <- tran$name()
  return(out)
}



