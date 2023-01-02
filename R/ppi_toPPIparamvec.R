# @title To or From vector form of parameters for PPI
# purely for testing hardcoded ppi estimators - returns canonical parameters
toPPIcannparam <- function(ALs, bL, beta, manifold = "sphere"){
  if (manifold == "sphere"){out <- ppi_paramvec(AL = ALs, bL = bL, beta = 1 + 2 * beta)}
  else if (manifold == "simplex"){out <- ppi_paramvec(AL = ALs, bL = bL, beta = beta)}
  else {stop("manifold not supported")}
  return(out)
}

# returns the length of the theta vector for given dimension p
ppithetalength <- function(p){
  p + #the beta
  (p-1) + #the diagonal of AL
  (p-2) * (p-1)/2 + #the upper triangle of AL
  (p-1) #the bL
}
ppiltheta2p <- function(ltheta){#ltheta is length of theta
  #ltheta = p + (p-1) + (p-2) * (p-1)/2 + (p-1)
  # = 3p - 2 + (p^2 - 3p + 2)/2
  # 0 = p^2/2 + 1.5p -1-ltheta
  # 0 = p^2 + 3p - (2+2ltheta)
  #p = (-3 pm sqrt(9 + 4(2+2ltheta)) / 2
  #p = (-3 + sqrt(9 + 8 + 8ltheta)) /2
  p <- (-3 + sqrt(9 + 8 + 8 * ltheta))/2
  return(p)
}

