# @title Internal function for building Score-Matching Objective Tapes
#' @param utape A measurment to use for taping
#' @param intheta A vector of parameters. NA values will be estimated, non-NA values will be fixed.
#' @param thetatape_creator A function that generates tape values for theta. Must take a single argument, `n` the number for values to generate
#' @param manifoldname Manifold with tranformation name ('Ralr', 'simplex'...)
#' @param llname Name of the log-likelihood function
#' @param utape An example observation (a single vector) to use for taping. The results shouldn't depend on `utape` so long as `utape` is in the manifold and the other functions are properly implemented.
#' @param usertheta A vector of parameter elements for the likelihood function. `NA` elements are to be fitted. Other elements are fixed at the provided value.
#' @param weightname The name of the divergence weight function ('ones' for manifolds without boundary).
#' @param acut The threshold `c` in the divergence weight function.
#' @param thetatape_creator A function that accepts an integer `n`, and returns a vector of `n` length. The function is used to fill in the `NA` elements of `usertheta` when building the tapes.
#' @param verbose If `TRUE` more details are printed when taping.
#' @description
#' Generates `CppAD` tapes for a model.
#' @details
#' When using, `CppAD` one first creates *tapes* of functions. These tapes can then be used for evaluating the function and its derivatives, and generating further tapes through argument swapping, differentiation and composition.
#' The taping relies on specifying typical argument values for the functions, so the programming is simplest when the function is defined without conditions.
#' Tapes can have both *independent* variables and *dynamic* variables.
#' The differentiation with `CppAD` occurs with respect to the independent variables.
#' Tapes of tapes are possible, including tapes that swap the independent and dynamic variables - this is how this package differentiates with respect to a dynamic variables.
#'
#' To build a tape for the score matching objective function, the package first tapes the map from a point \eqn{z} on the manifold to the value of the log-likelihood, where the independent variable is the \eqn{z}, and the dynamic variable is a vector of the *non*-fixed parameter values.
#' This tape is then used to generate a tape for the score matching objective with the non-fixed parameter values as the independent variable.
#' @return
#' Returns a list of the log-likelihood tape (note that the *input* for this tape is a function on the manifold), the tape of the score matching objective (the *input* here is the non-fixed parameter values), and some information used to generate the tape.
#' @examples
#' p <- 3
#' u <- movMF::rmovMF(1, rep(1, p))
#' ltheta <- p-1 + (p - 1) * p/2 + p #length of Fisher-Bingham parameter vector
#' intheta <- rep(NA, length.out = ltheta)
#' tapes <- buildsmotape("Snative", "FB",
#'               u, intheta,
#'               "ones", 1, verbose = FALSE
#'               )
#' pForward0(tapes$lltape, u, runif(n = ltheta))
#' pForward0(tapes$smotape, runif(n = ltheta), u)
# @export
buildsmotape <- function(manifoldname, llname,
                         utape, usertheta,
                         weightname = "ones", acut = 1,
                         thetatape_creator = function(n){seq(length.out = n)},
                         verbose = FALSE){
  starttheta <- t_u2s(usertheta, filler = thetatape_creator)
  isfixed <- t_u2i(usertheta)

  out <- buildsmotape_internal(manifoldname, llname,
                         utape, starttheta, isfixed,
                         weightname = weightname, acut = acut,
                         filler = thetatape_creator,
                         verbose = verbose)
  return(out)
}

buildsmotape_internal <- function(manifoldname, llname,
                         utape, starttheta, isfixed,
                         weightname = "ones", acut = 1,
                         filler = function(n){seq(length.out = n)},
                         verbose = FALSE){
  if(all(isfixed)){stop("All elements of theta are fixed")}
  thetatape <- starttheta

  if (!(manifoldname %in% c("simplex", "sphere"))){
    if (weightname != "ones"){warning("Manifold supplied has no boundary. Using weightname = 'ones' is strongly recommended.")}
  }
  if ((weightname == "ones") && (abs(acut - 1) > 1E-8)){
    warning("The value of 'acut' is ignored for weightname == 'ones'")
  }

  pman <- pmanifold(manifoldname)
  ztape <- ptoM(pman, utape) #the value of utape transformed to the manifold
  lltape <- ptapell(ztape, starttheta,
                    llname = llname, pman = pman,
                    fixedtheta = isfixed, verbose = verbose)
  stopifnot(is.numeric(acut))
  smotape <- ptapesmo(utape, t_si2f(starttheta, isfixed),
                      lltape, pman,
                      weightname, acut, verbose = verbose)
  return(list(
    lltape = lltape,
    smotape = smotape,
    info = list(
      name = llname,
      manifold = manifoldname,
      ulength = length(utape),
      starttheta = starttheta,
      isfixed = isfixed,
      weightname = weightname,
      acut = acut
    )
  ))
}
