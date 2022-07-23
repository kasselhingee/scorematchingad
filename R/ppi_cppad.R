
ppi_cppad <- function(prop, usertheta,
                      bdrythreshold = 1E-10, shiftsize = bdrythreshold, approxorder = 10,
                      pow = 1, man, weightname = hsqfun, acut = NULL, control = default_Rcgmin(), hsqfun = NULL){
  # process inputs
  stopifnot("matrix" %in% class(prop))
  p = ncol(prop)

  theta <- usertheta

  if (!(man %in% c("simplex", "sphere"))){
    if (weightname != "ones"){warning("Manifold supplied has no boundary. Setting weightname to 'ones'.")}
  }
  if (weightname == "ones"){
    if (!is.null(acut)){warning("The value of 'acut' is ignored for weightname == 'ones'")}
    acut <- 1 #set just for passing to Cpp
  }

  # prepare tapes
  tapes <- buildsmotape(man, "ppi",
                rep(1/p, p), theta,
                weightname = weightname,
                acut = acut, verbose = FALSE
                )

  # split data into boundary and interior
  datasplit <- simplex_boundarysplit(prop, bdrythreshold = bdrythreshold, shiftsize = shiftsize)

  opt <- cppadest(tapes$smotape, rep(0.2, sum(is.na(theta))), datasplit$interior,
               uboundary = datasplit$uboundary, boundaryapprox = datasplit$boundaryapprox,
               approxorder = approxorder,
               control = control)

  #process the theta and SE
  fixedtheta <- !is.na(theta)
  thetaest <- theta
  thetaest[!fixedtheta] <- opt$par
  SE <- fixedtheta * 0
  SE[!fixedtheta] <- opt$SE

  # make output
  list(
    prop = prop,
    est = c(list(theta = thetaest),
            fromPPIparamvec(thetaest, p)),
    SE = c(list(theta = SE),
           fromPPIparamvec(SE, p)),
    smval = opt$value,
    sqgradsize = opt$sqgradsize,
    counts = opt$counts,
    convergence = opt$convergence,
    message = opt$convergence
    )
}


