
ppi_cppad <- function(prop, stheta, isfixed,
                      bdrythreshold = 1E-10, shiftsize = bdrythreshold, approxorder = 10,
                      pow = 1, man, weightname = hsqfun, acut = NULL, control = default_Rcgmin(), hsqfun = NULL){
  # process inputs
  stopifnot("matrix" %in% class(prop))
  p = ncol(prop)

  if (!(man %in% c("simplex", "sphere"))){
    if (weightname != "ones"){warning("Manifold supplied has no boundary. Setting weightname to 'ones'.")}
  }
  if (weightname == "ones"){
    if (!is.null(acut)){warning("The value of 'acut' is ignored for weightname == 'ones'")}
    acut <- 1 #set just for passing to Cpp
  }

  # prepare tapes
  tapes <- buildsmotape_internal(man, "ppi",
                rep(1/p, p), 
                starttheta = stheta,
                isfixed = isfixed,
                weightname = weightname,
                acut = acut, verbose = FALSE
                )

  # split data into boundary and interior
  datasplit <- simplex_boundarysplit(prop, bdrythreshold = bdrythreshold, shiftsize = shiftsize)

  opt <- cppadest(tapes$smotape, t_si2f(stheta, isfixed), datasplit$interior,
               uboundary = datasplit$uboundary, boundaryapprox = datasplit$boundaryapprox,
               approxorder = approxorder,
               control = control)

  #process the theta and SE
  thetaest <- t_sfi2u(opt$par, stheta, isfixed)
  SE <- t_sfi2u(opt$SE, rep(0, length(stheta)), isfixed)

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


