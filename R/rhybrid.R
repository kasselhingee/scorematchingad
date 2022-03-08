#' @title Generate random observations from the PPI model
#' @description Given parameters of the PPI model, generates independent samples.
#' @param n Sample size
#' @param p Dimension (number of components)
#' @param beta0 The \eqn{\beta_0}{beta0} shape parameter vector
#' @param ALs The \eqn{A_L} parameter matrix
#' @param bL The \eqn{b_L} parameter vector
#' @param maxden This is the constant \eqn{log(C)} in (Scealy and Wood, 2021; Appendix A.1.1)
#' @return A list. The first element is the sample in the form of a matrix with `n` rows and `p` columns.
#' The second element is the maxden updated based on whether the sample exceeds the input maxden.
#' @examples
#' n=100
#' p=3
#' beta0=matrix(-0.8,p,1)
#' beta0[p]=-0.5
#'
#' muL=matrix(0,p-1,1)
#' muL[1:sum(p,-1)]=0.12
#' aa=matrix(1/500,p-1,1)
#' D=diag(as.vector(aa))
#' SigA=D
#' SigA[1,1]=SigA[1,1]*2
#' cor=0.5
#' SigA[1,2]=cor*sqrt(SigA[1,1]*SigA[2,2])
#' SigA[2,1]=SigA[1,2]
#' ALs=-0.5*solve(SigA)
#' bL=solve(SigA)%*%muL
#'
#' rhybrid(n,p,beta0,ALs,bL,4)
#' cdabyppi:::rhybrid_block(n,p,beta0,ALs,bL,4)
#' profvis::profvis({rhybrid(n,p,beta0,ALs,bL,4)})
#'
#' @export
rhybrid <- function(n,p,beta0,ALs,bL,maxden){
  # first simulate starting with a block of Dirichlet samples of size n.
  firstaccepted <- rhybrid_block(n,p,beta0,ALs,bL,maxden)
  maxden <- firstaccepted$maxden
  samples <- firstaccepted$accepted
  propaccepted <- max(nrow(samples) / n, 1E-3)
  # based on the number of samples that came out, simulate the remaining
  while (nrow(samples) < n){
    newsamples <- rhybrid_block(ceiling((n - nrow(samples)) * 1/propaccepted),p,beta0,ALs,bL,maxden)
    maxden <- newsamples$maxden
    samples <- rbind(samples, newsamples$accepted)
    # continue until n or more samples accepted
  }

  samples <- samples[1:n, ] # remove extra samples

  return(list(samp3=samples,maxden=maxden))
}


# simulates samples one at a time
rhybrid_singly <- function(n,p,beta0,ALs,bL,maxden)
{

	alpha=beta0+1
	coun=0
	samp1=matrix(0,1,p)
	count2=0
	tbL <- t(bL)
	while (coun < n)
	{

		count2=count2+1
		Uni=MCMCpack::rdirichlet(1, alpha)
		u=stats::runif(1,0,1)
		Uni_nop <- Uni[1:(p-1)]
		tUni_nop <- t(Uni[1:(p-1)])
		num=tUni_nop%*%ALs%*%Uni_nop + tbL%*%Uni_nop - maxden
		if (num > 0){maxden=num+maxden}
		#print(maxden)
		if (u < exp(num)){samp1=rbind(samp1,Uni);coun=coun+1}
		#print(coun)
	}
	samp3=samp1[2:length(samp1[,1]),]

	return(list(samp3=samp3,maxden=maxden))
}



# try generating samples without a 'while' loop
rhybrid_block <- function(n,p,beta0,ALs,bL,maxden){
  Uni <- MCMCpack::rdirichlet(n, beta0+1)
  Uni_nop <- Uni[, -p]
  nums <- .rowSums((Uni_nop %*% ALs) * Uni_nop, n, p-1) + #uT * ALs * u
    as.vector(Uni_nop %*% bL)  - maxden#t(bL) * u

  # update maxdens
  max_num <- max(nums)
  if (max_num > 0) {maxden=max_num+maxden}

  # accept some of points
  unif <- stats::runif(n,0,1)
  accepted <- Uni[unif < exp(nums), , drop = FALSE]

  return(list(accepted = accepted, maxden = maxden))
}
