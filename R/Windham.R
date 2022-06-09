#' @title Robust score matching estimates for the generalised-gamma form of the PPI model
#' @description Uses Windham weights after alr transform to estimate parameters for the PPI model with b_L=0.
#' The final (pth) component of beta0 is not estimated.
#' @param prop compositional data (n by p matrix)
#' @param cW the robustness tuning constant c
#' @param ALs_est initial values of A_L parameter matrix
#' @param bL_est value of b_L (this should be set to zero since not estimated)
#' @param beta0_est initial values of beta (beta[p]=0 and not estimated)
#' @param ind_weightA Roughly: the dimensions which have negative beta elements
#' @export
windham_diff=function(prop,cW,ALs_est,bL_est,beta0_est)
{
        stopifnot(all(bL_est == 0))
        p <- ncol(prop)
	sp=p-1
	x=c(1:sp)
	ind=combn(x, 2, FUN = NULL, simplify = TRUE)
	qind=length(ind[1,])

	stop1=0
	while(stop1==0)
	{
		#removing beta from weights
		d=-cW*beta0_est[1:sp]
                # the weights will be calculated directly from the ppi model density,
                # below modifies the A_L matrix to have zeros for the components that are not concentrated near zero
                ##### create the A_KK matrix from ALs_est. Elements of A_KK will be zero if dA is non-zero #####
		#removing A from weights
		dA=-cW*ALs_est
		for (j in 1:sp)
		{
			#putting some A's back in weights
			for (k in 1:sp)
			{
				if (ind_weightA[j]==0 && ind_weightA[k]==0){dA[j,k]=0} 
			}
		}
	
		weight_vec=matrix(0,n,1)
		
		ALs_estW=ALs_est
		for (j in 1:sp)
		{
			for (k in 1:sp)
			{
				if (dA[j,k]!=0){ALs_estW[j,k]=0}
			}
		}
                ###### A_KK creation finished #####
                # create the vector of weights from cW and ALs_estW (which is A_KK in Notes5.pdf)
                wwpar <- ppiparforww(beta0, ALs_est, bL_est, 1-ind_weightA)
                logden <- qldppi(prop, wwpar$beta0, wwpar$ALs, wwpar$bL)
                weight_mult=1
                weight_vec <- weight_mult * exp(cW*logden)
		weight_vec=n*(weight_vec/sum(weight_vec))

		#calculate scoring estimate:
		estimator=estimatorlog_weight(prop,0,weight_vec)$ppi
		estimate5=estimator

		previous1=ALs_est
		previous2=beta0_est

		#vectorise dA
		dA_vec=matrix(0,sum(sp,qind),1)
		for (j in 1:sp)
		{
			dA_vec[j]=dA[j,j]
		}
		for (j in 1:qind)
		{
			dA_vec[sum(j,sp)]=dA[ind[1,j],ind[2,j]]
			dA_vec[sum(j,sp)]=dA[ind[2,j],ind[1,j]]
		}
		as.vector(dA_vec)

		#update parameters
		estimate5[1:sum(sp,qind)]=(estimate5[1:sum(sp,qind)]-dA_vec)/(cW+1)
		estimate5[sum(sp,qind,1):sum(sp,qind,sp)]=(estimate5[sum(sp,qind,1):sum(sp,qind,sp)]-d)/(cW+1)
		x=c(1:sp)
		ind=combn(x, 2, FUN = NULL, simplify = TRUE)
		qind=length(ind[1,])
		ALs_est=matrix(0,sp,sp)
		for (j in 1:sp)
		{
			ALs_est[j,j]=estimate5[j]
		}
		for (j in 1:qind)
		{
			ALs_est[ind[1,j],ind[2,j]]=estimate5[sum(j,sp)]
			ALs_est[ind[2,j],ind[1,j]]=estimate5[sum(j,sp)]
		}
	
		beta0_est[1:sp]=estimate5[sum(sp,qind,1):sum(sp,qind,sp)]
	
		if ( (abs(beta0_est[1]-previous2[1])) < 0.000001 ) {stop1=1}


	}





	return(list(ALs_est=ALs_est,beta0_est=beta0_est,weight_vec=weight_vec,estimate5=estimate5))


}
		
ppiparforww <- function(beta0, ALs, bL, incinAL){
ALs_ww <- matrix(0, nrow(ALs), ncol(ALs))
ALs_ww[incinAL, incinAL] <- ALs[incinAL, incinAL]
return(list(beta0 = 0 * beta0, ALs = ALs_ww, bL = 0 * bL))
}
