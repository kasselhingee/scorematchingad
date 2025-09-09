#' @title Build a Tape of the Score Matching Integrand
#' @family tape builders
#' @family generic score matching tools
#' @param manifold A manifold object or name of a manifold used for estimation (which changes the score matching integrand)
#' @param uld Tape of an unnormalised log-density from [`tape_uld_inbuilt()`] or [`tape_uld()`]. For convenience, if `uld` is a string, it is passed to [`tape_uld_inbuilt()`].
#' @param transform If the `uld` is defined on a different embedding to `manifold`, then the name of the transformation that maps to the embedding of `manifold`.
#' @param xtape Vector to use for taping the `uld`. If there is a `transform`, then the transform is also applied to `xtape`. If omitted, the value of `uld$xtape` is used.
#' @param fixedparams A vector of fixed parameters for (re)taping `uld`. `NA` elements will remain *dynamic parameters*. Other elements will be fixed at the provided value. The same behaviour can also be obtained with [`fixdynamic()`] and [`fixindependent()`]. If not specified, no parameters of `uld` are fixed. For getting the correct dimensions, `fixedparams` is required if `uld` is a string, however there is no checking that `fixedparams` is the correct length for `uld` and `xtape`.
#' @param bdryw The tape of a boundary weight function from [`tape_bdryw()`] or [`tape_bdryw_inbuilt()`]. For convenience use `bdryw="ones"` for manifolds without boundary. If `bdryw="prodsq"` or `bdryw="minsq"` then boundary weight function for the simplex and positive orthant of the sphere are created with default threshold - see [`tape_bdrw_inbuilt()`] for more information.
#' @param dynparam_filler Optional. Defaults to values of `uld$dyntape` unless `uld` is a string. A function that generates values to be used for taping the dynamic parameters remaining after `fixedparams` is applied. Must take a single argument, `n` the number for values to generate. When `uld` is a string, then defaults to `function(n){seq(length.out = n)}`.
#' @param verbose If `TRUE` more details are printed when taping. These details are for debugging and will likely be comprehensible only to users familiar with the source code of this package.
#' @description
#' For a parametric model family, the function `tape_smi()` generates `CppAD` tapes for the unnormalised log-density of the model family and of the score matching integrand \eqn{A(z) + B(z) + C(z)} (defined in [`scorematchingtheory`]).
#' Three steps are performed by `tape_smi()`: first an object that specifies the manifold and any transformation from another manifold is created; then a tape of the unnormalised log-density is created; finally a tape of \eqn{A(z) + B(z) + C(z)} is created.
#' @details
#' If there is a `transform`, then `tape_smi()` first retapes via [`reembed()`] the `uld` for an independent variable \eqn{z} on `manifold` using the Jacobian provided in `transform` to account for the change in embedding (including the change of the uniform density/measure).
#' This tape has an independent variable \eqn{z} on `manifold`.
#' If supplied, `fixedparams` is applied using [`fixdynamic()`].
#' This tape is returned as `uld_reembed`.
#' 
#' This tape is then used to generate a tape for the score matching integrand where the parameters to estimate become the independent variable. 
#' The dynamic parameters of the tape correspond to the original manifold (domain) of `uld`, which reflects the fact that the score matching integrand is altered by `transform`, not the model or data (the data is not transformed and the model is unchanged).
#' This tape is returned as `smi`.
#' @references \insertAllCited{}

#' @return
#' A list of:
#'   + an [`Rcpp_ADFun`] object containing a tape of the unnormalised log-density using the metric of the "`end`" manifold (that is the independent variable is on the `end` manifold).
#'   + an [`Rcpp_ADFun`] object containing a tape of the score matching integrand with the non-fixed parameters of the model as the independent variable, and the measurements on the `end` manifold as the dynamic parameter.
#'   + some information about the tapes
#'
#' 
#' @examples
#' p <- 3
#' tapes <- tape_smi(manifold = "sph", 
#'                   uld = "vMF",
#'                   xtape = rep(1/sqrt(p), p),
#'                   fixedparams = rep(NA, p),
#'                   )
#' evaltape(tapes$uld_reembed, rep(1/sqrt(p), p), runif(n = p))
#' evaltape(tapes$smi, runif(n = p), rep(1/sqrt(p), p))
#' 
#' u <- rep(1/3, 3)
#' tapes <- tape_smi(manifold = "sph",
#'               uld = "ppi",
#'               transform = "sqrt"
#'               xtape = u,
#'               fixedparams = ppi_paramvec(p = 3), #no fixed params
#'               bdryw = "minsq" #uses default threshold for minsq
#'               )
#' evaltape(tapes$uld_reembed, sqrt(u), rppi_egmodel(1)$theta)
#' evaltape(tapes$smi, rppi_egmodel(1)$theta, u)
#' @export
tape_smi <- function(manifold,
                     uld,
                     transform = "none",
                     xtape = uld$xtape,
                     fixedparams = NA_real_*uld$dyntape,
                     bdryw = "ones",
                     dynparam_filler = NULL,
                     verbose = FALSE) {
  # if manifold a name, build manifold here
  if (is.character(manifold)){manifold <- make_manifold(manifold)}
  stopifnot(inherits(manifold, "Rcpp_man_ad"))

  # prepare uld
  if (all(!is.na(fixedparams))){stop("No elements of fixedparams are NA so no parameters will be estimated by this score matching.")}
  if (is.character(uld)){ #if uld is a name, build uld using xtape, fixedparams, and dynparam_filler
    stopifnot(is.numeric(xtape))
    stopifnot(length(xtape) > 0)
    stopifnot(length(fixedparams) > 0) #score matching doesnt make sense if no parameters to optimise
    if (is.null(dynparam_filler)){dynparam_filler <- function(n){seq(length.out = n)}}
    uld <- tape_uld_inbuilt(uld, 
              x = xtape, 
              theta = t_u2s(fixedparams, filler = dynparam_filler))
  } else { #else just fix parameters of passed uld using fiexdparams and dynparam_filler
    stopifnot(inherits("Rcpp_ADFun"))
    if (is.null(dynparam_filler)){
      theta <- fixedparams
      theta[t_u2i(fixedparams)] <- uld$dyntape[t_u2i(fixedparams)]
    } else {
      theta <- t_u2s(fixedparams, filler = dynparam_filler)
    }
    uld <- fixdynamic(uld, theta = theta, isfixed = t_u2i(usertheta)) # fix some of the model parameters if applicable
  }

  #reembed uld
  transform <- make_transform(transform)
  uld <- reembed(uld, tran = transform) #change the underlying metric of the manifold by using a different isometric embedding

  # choose between a canned boundary weight function or a custom boundary weight function
  if (is.character(bdryw)){ 
    bdryw <- tape_bdryw_inbuilt(bdryw, x = uld$xtape, acut = 0.2) #0.2 chosen as default for the 20% discussed around
  } else if (!inherits(bdryw, "Rcpp_ADFun")){
    stop("bdryw must be a name or a taped boundary weight function")
  }

  # combine everything for a full smi tape
  smitape <- tapesmd(uldtape = uld,
                        tran = transform,
                        M = manifold,
                        bdrywtape = bdryw,
                        verbose = verbose)
  return(list(
    uld_reembed = uld,
    smi = smitape,
    info = list(
      transform = transform$name(),
      manifold = manifold$name(),
      bdryw = bdryw
    )
  ))
}


