#' Make example models
#' @param n Number of samples to generate
#' @param maxden the log(C) maximum in simulating a PPI model
#' @param betaadd A number to add to the beta0 parameters - to changing the skewness of the marginal distributions
#' @return A list of the samples and model parameter information
#' @examples
#' sec2_3model(10, 4)
#' sec2_3model_p4(1000, 8)
#' @export
sec2_3model <- function(n, maxden = 4, betaadd = 0){
  #dimension
  p=3

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
  beta0 = beta0 + betaadd

  #simulate sample from PPI model
  samp1=cdabyppi:::rhybrid(n,p,beta0,ALs,bL, maxden)

  #simulated sample:
  samp3=samp1$samp3

  return(list(
    sample = samp3,
    p = p,
    beta0 = beta0,
    ALs = ALs,
    bL = bL,
    theta = toPPIparamvec(ALs, bL, beta0)
  ))
}

#' @describeIn sec2_3model
#' @export
sec2_3model_p4 <- function(n, maxden = 8){
  #dimension
  p=4

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
  theta <- c(diag(ALs), ALs[upper.tri(ALs)], bL, beta0)

  #simulate sample from PPI model
  samp1=cdabyppi:::rhybrid(n,p,beta0,ALs,bL,maxden)

  #simulated sample:
  samp3=samp1$samp3

  return(list(
    sample = samp3,
    p = p,
    beta0 = beta0,
    ALs = ALs,
    bL = bL,
    theta = toPPIparamvec(ALs, bL, beta0)
  ))
}
