#' @describeIn ppi_egmodel Simulates the 5-component model of \insertCite{@Section A.11, @scealy2022sc}{scorecompdir}. This model is a generalised gamma (i.e. \eqn{b_L} = 0).
#' The covariance-like matrix \eqn{A_L} is given in \insertCite{@Table 3ection A.11, @scealy2022sc}{scorecompdir}.
#' The parameters of the Dirichlet part are \eqn{\beta = (-0.80, -0.85, 0, -0.2, 0)}.
#' @references
#' \insertAllCited{}
#' @export
ppi_eggengamma <- function(n, maxden = 0){
  p <- 5
  ALs=matrix(0,p-1,p-1)
  bL=matrix(0,p-1,1)
  ALs[1,1]= -127480.0929
  ALs[1,2]= 14068.39057
  ALs[1,3]= 1782.261826
  ALs[1,4]=  -240.076568
  ALs[2,1]= ALs[1,2]
  ALs[2,2]= -8191.17253
  ALs[2,3]=  -8.002680
  ALs[2,4]= 374.693979
  ALs[3,1]= ALs[1,3]
  ALs[3,2]= ALs[2,3]
  ALs[3,3]= -46.638659
  ALs[3,4]= 9.027633
  ALs[4,1]= ALs[1,4]
  ALs[4,2]=  ALs[2,4]
  ALs[4,3]=  ALs[3,4]
  ALs[4,4]= -39.208915
  beta0=matrix(0,p,1)
  beta0[1]=-0.80
  beta0[2]=-0.85
  beta0[3]=0
  beta0[4]=-0.2
  beta0[5]=0

  #simulate sample from PPI model
  samp3=rppi(n,beta = beta0,AL = ALs,bL = bL, maxden=maxden)

  out <- list(
    sample = samp3,
    p = p,
    paramvec = toPPIparamvec(AL = ALs, bL = bL, beta = beta0),
    AL = ALs,
    bL = bL,
    beta = beta0)
  return(out)
}

