#' @title Generate random observations from the PPI model
#' @description Given parameters of the PPI model, generates independent samples.
#' @param n Sample size
#' @param p Dimension (number of components)
#' @param beta0 The \eqn{\beta_0}{beta0} shape parameter vector
#' @param ALs The \eqn{A_L} parameter matrix
#' @param bL The \eqn{b_L} parameter vector
#' @param maxden This is the constant \eqn{log(C)} in (Scealy and Wood, 2021; Appendix A.1.1)
#' @return A matrix with `n` rows and `p` columns and the maxden updated based on whether the sample exceeds the input maxden.
#' @export
rhybrid <- function(n,p,beta0,ALs,bL,maxden)
{

	alpha=beta0+1
	coun=0
	samp1=matrix(0,1,p)
	count2=0
	while (coun < n)
	{

		count2=count2+1
		Uni=MCMCpack::rdirichlet(1, alpha)
		u=stats::runif(1,0,1)
		num=t(Uni[1:sum(p,-1)])%*%ALs%*%t(t(Uni[1:sum(p,-1)]))+t(bL)%*%t(t(Uni[1:sum(p-1)]))-maxden
		if (num > 0){maxden=num+maxden}
		#print(maxden)
		if (u < exp(num)){samp1=rbind(samp1,Uni);coun=coun+1}
		#print(coun)
	}
	samp3=samp1[2:length(samp1[,1]),]

	return(list(samp3=samp3,maxden=maxden))
}

