#' @title Build and simulate from example PPI models
#' @param n Number of samples to generate
#' @param maxden the log(C) maximum in simulating a PPI model
#' @param betaadd A number to add to the \eqn{\beta} parameter vector, which may be useful for experimenting with how 'Dirichlet' the model look.
#' @return A list of the samples and model parameter information
#' @examples
#' ppi_egmodel(1000)
#' @description 
#' The model is the 3-component model from \insertCite{@Section 2.3, @scealy2022sc}{cdabyppi}.
#' @references
#' \insertAllCited{}
#' @export
ppi_egmodel <- function(n, maxden = 4, betaadd = 0){
  mats <- pars_sec2dot3model(3)
  mats$beta = mats$beta + betaadd

  #simulate sample from PPI model
  samp1=rppi(n,3,beta0 = mats$beta,ALs = mats$AL,bL = mats$bL, maxden)

  #simulated sample:
  samp3=samp1$samp3

  out <- c(list(
    sample = samp3,
    p = 3,
    theta = toPPIparamvec(mats$AL, mats$bL, mats$beta)
    ),
    mats,
    list(beta0 = mats$beta))


  return(out)
}

# @describeIn ppi_egmodel
# @export
ppi_egmodel_p4 <- function(n, maxden = 8){
  mats <- pars_sec2dot3model(4)

  #simulate sample from PPI model
  samp1=rppi(n,4,beta0 = mats$beta,ALs = mats$AL,bL = mats$bL, maxden)

  #simulated sample:
  samp3=samp1$samp3
  
  out <- c(list(
    sample = samp3,
    p = 3,
    theta = toPPIparamvec(mats$AL, mats$bL, mats$beta)
    ),
    mats,
    list(beta0 = mats$beta))


  return(out)
}


# internal function for building Section 2.3-like model for any number of components
# it pretty poor for anything but p = 3
# returns the parameter matrices/vectors
pars_sec2dot3model <- function(p){
  #parameters for the PPI model
  muL=matrix(0,p-1,1)
  muL[1:sum(p,-1)]=0.12
  aa=matrix(1/500,p-1,1)
  D=diag(as.vector(aa))
  SigA=D
  SigA[1,1]=SigA[1,1]*2
  cor=0.5
  SigA[1,2]=cor*sqrt(SigA[1,1]*SigA[2,2])
  SigA[2,1]=SigA[1,2]
  ALs=-0.5*solve(SigA)
  bL=solve(SigA)%*%muL
  beta0=matrix(-0.8,p,1)
  beta0[p]=-0.5
  return(list(
    AL = ALs,
    bL = bL,
    beta = beta0
  ))
}
