#' @title Internal function for building Score-Matching Objective Tapes
#' @examples
#' @param utape A measurment to use for taping
#' @param intheta A vector of parameters. NA values will be estimated, non-NA values will be fixed.
#' @param thetatape_creator A function that generates tape values for theta. Must take a single argument, `n` the number for values to generate
#' @examples
#' p <- 3
#' u <- Directional::rvmf(1, rep(1, p), 0) #uniform distribution on the sphere p-space
#' ltheta <- p-1 + (p - 1) * p/2 + p
#' intheta <- rep(NA, length.out = ltheta)
#' tapes <- buildsmotape("Snative", "FB",
#'               u, intheta,
#'               "ones", 1, verbose = FALSE
#'               )
#' pForward0(tapes$lltape, u, runif(n = ltheta))
#' pForward0(tapes$smotape, runif(n = ltheta), u)
#' @export
buildsmotape <- function(manifoldname, llname,
                         utape, intheta,
                         weightname = "ones", acut = 1,
                         thetatape_creator = function(n){seq(length.out = n)},
                         verbose = FALSE){
  fixedtheta <- !is.na(intheta)
  thetatape <- intheta
  thetatape[!fixedtheta] <- thetatape_creator(sum(!fixedtheta))

  if (!(manifoldname %in% c("simplex", "sphere"))){
    if (weightname != "ones"){warning("Manifold supplied has no boundary. Using weightname = 'ones' is strongly recommended.")}
  }
  if ((weightname == "ones") && (abs(acut - 1) > 1E-8)){
    warning("The value of 'acut' is ignored for weightname == 'ones'")
  }

  pman <- pmanifold(manifoldname)
  ztape <- ptoM(pman, utape) #the value of utape transformed to the manifold
  lltape <- ptapell(ztape, thetatape,
                    llname = llname, pman = pman,
                    fixedtheta = fixedtheta, verbose = verbose)
  stopifnot(is.numeric(acut))
  smotape <- ptapesmo(utape, thetatape[!fixedtheta],
                      lltape, pman,
                      weightname, acut, verbose = verbose)
  return(list(
    lltape = lltape,
    smotape = smotape,
    info = list(
      name = llname,
      manifold = manifoldname,
      ulength = length(utape),
      intheta = intheta,
      tapedtheta = thetatape,
      weightname = weightname,
      acut = acut
    )
  ))
}
