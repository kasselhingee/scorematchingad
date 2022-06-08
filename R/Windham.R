#######################################################
##estimatorlog_weight calculates the score matching estimates 
##using the log-ratio transformation
##for the generalised gamma model (A_L and beta estimated).
##note that beta[p] is fixed at a chosen value
##betap: chosen value of beta[p]
##prop: compositional data (n by p matrix)
##weightW: weight vector
#######################################################


estimatorlog_weight <- function(prop,betap,weightW)
{

	sp=p-1

	x=c(1:sp)
	ind=combn(x, 2, FUN = NULL, simplify = TRUE)
	qind=length(ind[1,])

	V1S_sum=matrix(0,sum(sp,sp,qind),1)
	V2S_sum=matrix(0,sum(sp,sp,qind),1)
	V3S_sum=matrix(0,sum(sp,sp,qind),sum(sp,sp,qind))


	for (b in 1:n)
	{

		V1A=matrix(0,sp,sp)
		for (j in 1:sp)
		{	
			for (i in 1:sp)
			{
				V1A[i,j]=-2*(prop[b,i]^2)*prop[b,j]+6*(prop[b,i]^2)*(prop[b,j]^2)
				V1A[j,i]=-2*(prop[b,j]^2)*prop[b,i]+6*(prop[b,j]^2)*(prop[b,i]^2)
			}
			V1A[j,j]=4*prop[b,j]^2-10*prop[b,j]^3+6*prop[b,j]^4
		}


		V1B=matrix(0,qind,sp)
		for (j in 1:qind)
		{

			for (k in 1:sp)
			{
				if (ind[1,j]==k){V1B[j,k]=2*prop[b,ind[1,j]]*prop[b,ind[2,j]]-12*(prop[b,ind[1,j]]^2)*prop[b,ind[2,j]]+12*(prop[b,ind[1,j]]^3)*prop[b,ind[2,j]]}
				else if (ind[2,j]==k){V1B[j,k]=2*prop[b,ind[2,j]]*prop[b,ind[1,j]]-12*(prop[b,ind[2,j]]^2)*prop[b,ind[1,j]]+12*(prop[b,ind[2,j]]^3)*prop[b,ind[1,j]]}
				else {V1B[j,k]=-4*prop[b,k]*(1-prop[b,k])*prop[b,ind[1,j]]*prop[b,ind[2,j]]+8*(prop[b,k]^2)*prop[b,ind[1,j]]*prop[b,ind[2,j]]}

			}

		}


		V1C=matrix(0,sp,sp)
		for (j in 1:sp)
		{

			V1C[,j]=-prop[b,j]*(1-prop[b,j])
		}



		V2A=matrix(0,sp,sp)
		for (j in 1:sp)
		{
	
			for (i in 1:sp)
			{
				V2A[i,j]=-2*(prop[b,i]^2)*prop[b,j]
				V2A[j,i]=-2*(prop[b,j]^2)*prop[b,i]
			}
			V2A[j,j]=2*(prop[b,j]^2)*(1-prop[b,j])
		}


		V2B=matrix(0,qind,sp)
		for (j in 1:qind)
		{

			for (k in 1:sp)
			{
				if (ind[1,j]==k){V2B[j,k]=2*prop[b,ind[1,j]]*prop[b,ind[2,j]]*(1-2*prop[b,ind[1,j]])}
				else if (ind[2,j]==k){V2B[j,k]=2*prop[b,ind[2,j]]*prop[b,ind[1,j]]*(1-2*prop[b,ind[2,j]])}
				else {V2B[j,k]=-4*prop[b,k]*prop[b,ind[1,j]]*prop[b,ind[2,j]]}

			}

		}


		V2C=matrix(0,sp,sp)
		for (j in 1:sp)
		{
	
			for (i in 1:sp)
			{
				V2C[i,j]=-1*prop[b,j]
				V2C[j,i]=-1*prop[b,i]
			}
			V2C[j,j]=(1-prop[b,j])
		}


		V1=rbind(V1A,V1B,V1C)
		V2=rbind(V2A,V2B,V2C)

		V1S=matrix(0,sum(sp,sp,qind),1)
		V2S=matrix(0,sum(sp,sp,qind),1)
		V3S=matrix(0,sum(sp,sp,qind),sum(sp,sp,qind))

		for (j in 1:sp)
		{	
			V1S=V1S+V1[,j]
			V2S=V2S+V2[,j]*prop[b,j]
			V3S=V3S+t(t(V2[,j]))%*%t(V2[,j])
		}

		V1S_sum=V1S_sum+V1S*weightW[b]
		V2S_sum=V2S_sum+V2S*weightW[b]
		V3S_sum=V3S_sum+V3S*weightW[b]


	}

	V1S_sum=V1S_sum/n
	V2S_sum=V2S_sum/n
	V3S_sum=V3S_sum/n

	alpha_p=betap+1

	d=((alpha_p*V2S_sum)-V1S_sum)
	W=V3S_sum
	ppi=solve(W)%*%d

	tot=sum(sp,qind,sp)
	ppi[sum(sp,qind,1):tot]=ppi[sum(sp,qind,1):tot]-1

	#ppi is the score matching estimate
	return(list(ppi=ppi,W=W,d=d))


}



#####################################################
##windham_diff calculates the robust score matching estimates
##using an iterative re-weighting algorithm. 
##prop: compositional data (n by p matrix)
##cW: the robustness tuning constant c
##ALs_est: initial values of A_L parameter matrix
##bL_est: value of b_L (this should be set to zero since not estimated)
##beta0_est: initial values of beta (beta[p]=0 and not estimated)
#####################################################


windham_diff=function(prop,cW,ALs_est,bL_est,beta0_est)
{

	sp=p-1
	x=c(1:sp)
	ind=combn(x, 2, FUN = NULL, simplify = TRUE)
	qind=length(ind[1,])

	stop1=0
	while(stop1==0)
	{
		#removing beta from weights
		d=-cW*beta0_est[1:sp]
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
		for (i in 1:n)
		{
			weight_mult=1
			weightA=weight_mult*exp(cW*(t(prop[i,1:sp])%*%ALs_estW%*%t(t(prop[i,1:sp])))+cW*(t(bL_est)%*%t(t(prop[i,1:sp]))))
			weight_vec[i]=weightA
		}
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
