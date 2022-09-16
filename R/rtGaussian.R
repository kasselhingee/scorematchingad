# @title Generate random observations from a truncated Gaussian model
# @description Generates random observations from a truncated Gaussian model with \eqn{\beta=0}{beta=0}.
# Purely for testing Model 4 from Section A.11 of Scealy and Wood (2022).
# @param n Sample size
# @param p Dimension (number of components)
# @param muL Untruncated Gaussian mean. See (Scealy and Wood, 2022; Appendix A.1.1).
# @param SigA Untruncated Gaussian covariance matrix. See (Sealy and Wood, 2022; Appendix A.1.1).
# @return A matrix with `n` rows and `p` columns.
# @export
rtGaussian <- function(n,p,muL,SigA)
{

	coun=0
	samp2=matrix(0,1,p-1)
	count2=0
	while (coun < n)
	{
		count2=count2+1
		rand=MASS::mvrnorm(1,mu=muL,Sigma=SigA)
		if (min(rand) >= 0 && sum(rand) < 1){samp2=rbind(samp2,rand);coun=coun+1}
		#print(coun)
	}
	samp2=samp2[2:length(samp2[,1]),]
	com=matrix(0,n,1)
	samp2=cbind(samp2,com)
	samp2[,p]=1-rowSums(samp2[,1:sum(p,-1)])
	samp3=samp2

	return(samp3)
}


