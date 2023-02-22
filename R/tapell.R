#' @title Tape a log-likelihood
#' @param pmanifoldtransform An `Rcpp::XPtr` object pointing to a manifold and transformation to the manifold from the natural space of the log-likelihood. Generate `pmanifoldtransform` with [`manifoldtransform()`].
#' @param llname The name of the log-likelihood function to tape
#' @param ytape An example measurement value to use for creating the tape. In the natural manifold of the log-likelihood function. `ytape` will be converted to the manifold according to `manifoldname` before taping. 
#' Please ensure that `ytape` is the interior of the manifold, and it is probably best if all components of `ytape` are non-zero when transformed form onto the manifold provided in `manifoldname`.
#' @param usertheta A vector of parameter elements for the likelihood function. `NA` elements will become *dynamic parameters*. Other elements will be fixed at the provided value. The length of `usertheta` must be the correct length for the log-likelihood, no checking is conducted.
#' @param thetatape_creator A function that accepts an integer `n`, and returns a vector of `n` length. The function is used to fill in the `NA` elements of `usertheta` when building the tapes. Please ensure that the values filled by `thetatape_creator` lead to plausible parameter vectors for the chosen log-likelihood.
#' @param verbose If `TRUE`, some information about the tape is printed.
#' @description Creates a `CppAD` tape of the log-likelihood as a function of values on the provided manifold. The log-likelihood accounts for the change of measure between the natural manifold of the log-likelihood and the manifold with transformation provided by `manifoldname`.
#' @return 
#' An `Rcpp::XPtr` object that points to a `CppAD::ADFun` object.
#' The object will also have attributes:
#'  + `fname` A name of the taped function, derived from `llname`
#'  + `ytape` The value of `ytape`
#'  + `dyntape` The value of non-fixed elements of `usertheta` used for taping.
#'  + `usertheta` The vector `usertheta`.
#' @examples 
#' simpwsqrt <- manifoldtransform("sphere")
#' ppitape <- tapell(llname = "ppi",
#'                   ytape = c(0.2, 0.3, 0.5),
#'                   usertheta = ppi_paramvec(p = 3), 
#'                   pmanifoldtransform = simpwsqrt)
#' pForward0(ppitape, 
#'   sqrt(rep(1/3, 3)), 
#'   ppi_paramvec(p = 3, AL=0, bL=0, beta=c(0,0,0.5)))
#' pJacobian(ppitape, 
#'   sqrt(rep(1/3, 3)), 
#'   ppi_paramvec(p = 3, AL=0, bL=0, beta=c(0,0,0.5)))
#' pHessian(ppitape, 
#'   sqrt(rep(1/3, 3)), 
#'   ppi_paramvec(p = 3, AL=0, bL=0, beta=c(0,0,0.5)))
#' @section Warning: There is limited checking of the inputs.
#' @export
tapell <- function(llname,
                   ytape,
                   usertheta,
                   pmanifoldtransform,
                   thetatape_creator = function(n){seq(length.out = n)},
                   verbose = FALSE){
  starttheta <- t_u2s(usertheta, filler = thetatape_creator)
  ztape <- ptoM(pmanifoldtransform, ytape) #the value of ytape transformed to the manifold
  lltape <- ptapell(ztape, starttheta,
                    llname = llname, 
                    pman = pmanifoldtransform,
                    fixedtheta = t_u2i(usertheta),
                    verbose = verbose)
  out <- ADFun$new(ptr = lltape,
                   name = paste(attr(pmanifoldtransform, "name"), llname, sep = "-"),
                   xtape = ztape,
                   dyntape =  starttheta[!t_u2i(usertheta)])

  attr(out, "ytape") <- ytape
  attr(out, "manifoldtransform") <- attr(pmanifoldtransform, "name")
  attr(out, "usertheta") <- usertheta
  return(out)
}



