# @title To or From vector form of parameters for PPI
# @export
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

combparam2uppertriorder <- function(ltheta){
  p <- ppiltheta2p(ltheta)
  #order of AL part
  inthetaorder <- combparam2uppertriorder_matrix(p-1 + (p-1)*(p-2)/2)
  # add rest of theta
  inthetaorder <- c(inthetaorder,
                    seq.int(max(inthetaorder) + 1, length.out = p - 1 + p))
  stopifnot(length(inthetaorder) == ltheta) #check that length produced is correct
  stopifnot(all(inthetaorder[inthetaorder] == seq(1, length.out = ltheta))) #check that ordering index is correct
  return(inthetaorder)
}


combparam2uppertriparam <- function(theta){
  inthetaorder <- combparam2uppertriorder(length(theta))
  return(theta[inthetaorder])
}
