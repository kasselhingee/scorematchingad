#' @title Load and simulate example PPI models
#' @param n Number of samples to generate
#' @param maxden the log(C) maximum in simulating a PPI model. See [`rppi()`].
#' @return A list: 
#'  * `sample` A matrix of the simulated samples (`n` rows)
#'  * `p` The number of components of the model
#'  * `theta` The PPI parameter vector
#'  * `AL` The \eqn{A_L} parameter matrix
#'  * `bL` The \eqn{b_L} parameter vector
#'  * `beta` The \eqn{\beta} parameter vector 
#' @examples
#' ppi_egmodel(1000)
#' ppi_eggengamma(100)
#' @description 
#' The function `ppi_egmodel()` simulates the 3-component PPI model from \insertCite{@Section 2.3, @scealy2022sc}{scorecompdir}. The function `ppi_eggengamma()` simulates from the 5-component model of \insertCite{@Section A.11, @scealy2022sc}{scorecompdir}, which is a generalised gamma model (i.e. \eqn{b_L} = 0).
#' Model parameters and the simulated sample are returned.
#' @references
#' \insertAllCited{}
#' @export
ppi_egmodel <- function(n, maxden = 4){
  mats <- pars_sec2dot3model(3)
  mats$beta = mats$beta

  #simulate sample from PPI model
  samp3=rppi(n,beta = mats$beta,AL = mats$AL,bL = mats$bL, maxden = maxden)

  out <- c(list(
    sample = samp3,
    p = 3,
    theta = ppi_paramvec(AL = mats$AL, bL = mats$bL, beta = mats$beta)
    ),
    mats,
    list(beta0 = mats$beta,
         ALs = mats$AL))


  return(out)
}

# @rdname ppi_egmodel
# @export
ppi_egmodel_p4 <- function(n, maxden = 8){
  mats <- pars_sec2dot3model(4)

  #simulate sample from PPI model
  samp3=rppi(n,beta = mats$beta,AL = mats$AL,bL = mats$bL, maxden = maxden)

  out <- c(list(
    sample = samp3,
    p = 4,
    theta = ppi_paramvec(AL = mats$AL, bL = mats$bL, beta = mats$beta)
    ),
    mats,
    list(beta0 = mats$beta,
         ALs = mats$AL))


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
