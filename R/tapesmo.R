#' @title Tape a score matching objective function
#' @description Generates a `CppAD` tape of the empirical score matching objective function for a single measurement (\eqn{\tilde\psi} in Score Matching vignette). Requires a tape of the log-likelihood function and the corresponding manifold with the transformation. 
#' @param pmanifoldtransform An `Rcpp::XPtr` object pointing to a manifold and transformation to the manifold from the natural space of the log-likelihood. Generate `pmanifoldtransformation` with [`pmanifold()`].
#' @param lltape Tape of the log-likelihood function from the manifold represented by `pmanifoldtransform`. Construct `lltape` using [`tapell()`].
#' @param divweight The name of the divergence weight function ("ones" for manifolds without boundary). For the simplex and positive orthant of the sphere, "prodsq" and "minsq" are possible.
#' @param acut The threshold \eqn{a_c} in the divergence weight function `divweight`. Ignored for `divweight = "ones"`.
#' @examples 
#' sqrtman <- pmanifold("sphere")
#' ppitape <- tapell(llname = "ppi",
#'                   xtape = c(0.2, 0.3, 0.5),
#'                   usertheta = ppi_paramvec(p = 3), 
#'                   pmanifoldtransform = sqrtman)
#' ppismotape <- tapesmo(lltape = ppitape,
#'                       pmanifoldtransform = sqrtman,
#'                       divweight = "minsq",
#'                       acut = 0.1,
#'                       verbose = FALSE)
#' pForward0(ppismotape, 
#'   ppi_paramvec(p = 3, AL=0, bL=0, beta=c(-0.1,-0.1,0.5)),
#'   rep(1/3, 3)) 
# objs <- buildsmotape("sphere", "ppi", c(0.2, 0.3, 0.5),
#                      weightname = "minsq", acut = 0.1,
#                      ppi_paramvec(p = 3))
# pForward0(objs$smotape, 
#   ppi_paramvec(p = 3, AL=0, bL=0, beta=c(-0.1,-0.1,0.5)),
#   c(0.1, 0.1, 0.8)) 
#' @export
tapesmo <- function(lltape,
                   pmanifoldtransform,
                   divweight,
                   acut = 1,
                   verbose = FALSE){
  stopifnot(is.numeric(acut))
  
  smotape <- ptapesmo(attr(lltape, "xtape"),
                      attr(lltape, "dyntape"),
                      lltape,
                      pmanifoldtransform,
                      divweight, 
                      acut, 
                      verbose = verbose)
  attr(smotape, "fname") <- paste(attr(lltape, "fname"), "smo", sep = "-")
  attr(smotape, "xtape") <- attr(lltape, "dyntape")
  attr(smotape, "usertheta") <- attr(lltape, "usertheta")
  attr(smotape, "dyntape") <- attr(lltape, "xtape")
  attr(smotape, "divweight") <- divweight
  attr(smotape, "acut") <- acut
  return(smotape)
}
