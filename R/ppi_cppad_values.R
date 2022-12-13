# OBSOLETE - JUST USED FOR TESTING
# A function that gives the sm val, gradient and hessian for a given theta
# arguments are the same as ppi_cppad()
# stheta is the value that is used for evaluating, this value is also used for taping the cppad functions
ppi_cppad_values <- function(prop, stheta, isfixed,
                      bdrythreshold = 1E-10, shiftsize = bdrythreshold, approxorder = 10,
                      pow = 1, man, weightname = hsqfun, acut = NULL, control = default_Rcgmin(), hsqfun = NULL,
                      w = NULL){

  p <- ncol(prop)

  # prepare tapes
  tapes <- buildsmotape_internal(man, "ppi",
                rep(1/p, p),
                starttheta = stheta,
                isfixed = isfixed,
                weightname = weightname,
                acut = acut, verbose = FALSE
                )

  # split data into boundary and interior
  datasplit <- simplex_boundarysplit(prop, bdrythreshold = bdrythreshold, shiftsize = shiftsize, w = w)

  numboundarypoints <- switch(1 + is.null(datasplit$uboundary), 
                           nrow(datasplit$uboundary),
                           0)

  eginteriorpt <- rep(1/p, p)
  theta_for_taping <- t_si2f(stheta, tapes$info$isfixed)
  Jsmofun <- pTapeJacobian(smotape, theta_for_taping, eginteriorpt)
  Hsmofun <- pTapeJacobian(Jsmofun, theta_for_taping, eginteriorpt)
  smofun_u <- swapDynamic(smotape, eginteriorpt, theta_for_taping)
  Jsmofun_u <- swapDynamic(Jsmofun, eginteriorpt, theta_for_taping)
  Hsmofun_u <- swapDynamic(Hsmofun, eginteriorpt, theta_for_taping)


  objval <- smobj(smofun = tapes$smotape,
        theta = t_si2f(stheta, tapes$info$isfixed),
        utabl = datasplit$interior,
        smofun_u = smofun_u,
        uboundary = datasplit$uboundary,
        boundaryapprox = datasplit$boundaryapprox,
        approxorder = approxorder,
        w = datasplit$winterior,
        wboundary = datasplit$wboundary
        )

  gradval <- smobjgrad(smofun = tapes$smotape,
        theta = t_si2f(stheta, tapes$info$isfixed),
        utabl = datasplit$interior,
        Jsmofun_u = Jsmofun_u,
        uboundary = datasplit$uboundary,
        boundaryapprox = datasplit$boundaryapprox,
        approxorder = approxorder,
        w = datasplit$winterior,
        wboundary = datasplit$wboundary
        )
  
  hessval <- smobjhess(smofun = tapes$smotape,
        theta = t_si2f(stheta, tapes$info$isfixed),
        utabl = datasplit$interior,
        Hsmofun_u = Hsmofun_u,
        uboundary = datasplit$uboundary,
        boundaryapprox = datasplit$boundaryapprox,
        approxorder = approxorder,
        w = datasplit$winterior,
        wboundary = datasplit$wboundary
        )
  
  return(list(obj = objval,
              grad = gradval,
              hess = hessval
              ))
}

