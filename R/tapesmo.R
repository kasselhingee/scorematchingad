#' @describeIn buildsmotape Generates a `CppAD` tape of the empirical score matching objective function for a single measurement (\eqn{\tilde\psi} in Score Matching vignette). Requires a tape of the log-likelihood function and the corresponding manifold with the transformation. 
#' @param tran A transform object (of type `Rcpp_transform_ad`), typically created by [`manifoldtransform()`].
#' @param man A transform object (of type `Rcpp_man_ad`) matching `tran`, and typically created by [`manifoldtransform()`].
#' @param lltape Tape of the log-likelihood function constructed using [`tapell()`].
#' @param divweight The name of the boundary weight function ("ones" for manifolds without boundary). For the simplex and positive orthant of the sphere, "prodsq" and "minsq" are possible.
#' @param acut The threshold \eqn{a_c} in the boundary weight function `divweight`. Ignored for `divweight = "ones"`.
#' @param verbose If `TRUE`, some information about the tape is printed.
#' @return 
#' `tapesmo()` returns an [`ADFun`] object.  
#' @examples 
#' sqrtman <- manifoldtransform("sim", "sqrt", "sph")
#' ppitape <- tapell(llname = "ppi",
#'                   ytape = c(0.2, 0.3, 0.5),
#'                   usertheta = ppi_paramvec(p = 3), 
#'                   tran = sqrtman$tran)
#' ppismotape <- tapesmo(lltape = ppitape,
#'                       tran = sqrtman$tran,
#'                       man = sqrtman$man,
#'                       divweight = "minsq",
#'                       acut = 0.1,
#'                       verbose = FALSE)
#' pForward0(ppismotape$ptr, 
#'   ppi_paramvec(p = 3, AL=0, bL=0, beta=c(-0.1,-0.1,0.5)),
#'   rep(1/3, 3)) 
# objs <- buildsmotape("sphere", "ppi", c(0.2, 0.3, 0.5),
#                      divweight = "minsq", acut = 0.1,
#                      ppi_paramvec(p = 3))
# pForward0(objs$smotape, 
#   ppi_paramvec(p = 3, AL=0, bL=0, beta=c(-0.1,-0.1,0.5)),
#   c(0.1, 0.1, 0.8)) 
#' @export
tapesmo <- function(lltape,
                   tran,
                   man,
                   divweight,
                   acut = 1,
                   verbose = FALSE){
  inherits(tran, "Rcpp_transform_ad")
  inherits(man, "Rcpp_man_ad")
  stopifnot(is.numeric(acut))
  
  smotape <- ptapesmo(attr(lltape, "ytape"),
                      lltape$dyntape,
                      lltape$ptr,
                      tran,
                      man,
                      divweight, 
                      acut, 
                      verbose = verbose)
  out <- ADFun$new(ptr = smotape,
                   name = paste(lltape$name, "smo", sep = "-"),
                   xtape = as.numeric(lltape$dyntape),
                   dyntape = as.numeric(attr(lltape, "ytape")),
                   usertheta = as.numeric(lltape$usertheta))
  return(out)
}
