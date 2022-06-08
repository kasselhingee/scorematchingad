######################################################################
##rhybrid generates random observations from the PPI model
##n: sample size
##p: dimension
##beta0: beta shape parameter vector
##ALs: A_L parameter matrix
##bL: b_L parameter vector 
##maxden: is the constant log(C) in Appendix A.1.3. in Scealy and Wood (2022)
######################################################################


library(coda)
library(MCMCpack)

rhybrid <- function(n,p,beta0,ALs,bL,maxden)
{

	alpha=beta0+1
	coun=0
	samp1=matrix(0,1,p)
	count2=0
	while (coun < n)
	{

		count2=count2+1
		Uni=rdirichlet(1, alpha)
		u=runif(1,0,1)
		num=t(Uni[1:sum(p,-1)])%*%ALs%*%t(t(Uni[1:sum(p,-1)]))+t(bL)%*%t(t(Uni[1:sum(p-1)]))-maxden
		if (num > 0){maxden=num+maxden}
		#print(maxden)
		if (u < exp(num)){samp1=rbind(samp1,Uni);coun=coun+1}
		#print(coun)
	}
	samp3=samp1[2:length(samp1[,1]),]

	return(list(samp3=samp3,maxden=maxden))
}



#######################################################
##estimator1 calculates Score1ac and Score1 estimators in 
##Scealy and Wood (2022)
##for the model with A_L and b_L estimated and beta fixed.
##prop: compositional data (n by p matrix)
##acut: a_c for h function (if acut > 1 then this is Score1)
##incb: if incb=1 then b_L is estimated otherwise b_L is fixed at zero (omitted)
#######################################################



estimator1 <- function(prop,acut,incb)
{

	#response on sphere scale
	z=sqrt(prop)

	#h function
	h=matrix(1,n,1)
	for (j in 1:p)
	{
		h=h*z[,j]
	}
	indh=matrix(0,n,1)
	h2=h
	for (j in 1:n)
	{
		indh[j]=1
		zmin=z[j,1]
	
		for (i in 2:p)
		{
			zmin_prev=zmin
			zmin=min(zmin,z[j,i])
			if (zmin_prev > zmin){indh[j]=i}
		}
		zmin_prev=zmin
		zmin=min(zmin,acut)
		if (zmin_prev > zmin){indh[j]=0}
		h2[j]=zmin
	
	}
	h=h2


	###################
	##calculate W
	##################

	sp=p-1

	A1=matrix(0,1,p-1)
	for (j in 1:sp)
	{
		A1[j]=mean(16*(h^2)*z[,j]^6)
	}

	x=c(1:sp)
	ind=combn(x, 2, FUN = NULL, simplify = TRUE)
	qind=length(ind[1,])


	C1=matrix(0,1,p-1)
	for (j in 1:sp)
	{
		C1[j]=mean(8*(h^2)*z[,j]^4)
	}

	AA=diag(as.vector(A1))


	AB=matrix(0,sp,qind)
	for (i in 1:sp)
	{
		for (j in 1:qind)
		{
			if (i==ind[1,j]){AB[i,j]=mean(16*z[,ind[1,j]]^4*z[,ind[2,j]]^2*h^2)}
			if (i==ind[2,j]){AB[i,j]=mean(16*z[,ind[1,j]]^2*z[,ind[2,j]]^4*h^2)}
		}
	}

	BB=matrix(0,qind,qind)
	for (i in 1:qind)
	{
		for (j in 1:qind)
		{
			if (ind[1,i]==ind[1,j] & ind[2,i]==ind[2,j]){BB[i,j]=mean(16*z[,ind[1,j]]^4*z[,ind[2,j]]^2*h^2)+mean(16*z[,ind[1,j]]^2*z[,ind[2,j]]^4*h^2)}
			else if (ind[1,i]==ind[1,j]){BB[i,j]=mean((h^2)*16*z[,ind[2,j]]^2*z[,ind[2,i]]^2*z[,ind[1,j]]^2)}
			else if (ind[2,i]==ind[2,j]){BB[i,j]=mean((h^2)*16*z[,ind[1,j]]^2*z[,ind[1,i]]^2*z[,ind[2,j]]^2)}
			else if (ind[1,i]==ind[2,j]){BB[i,j]=mean((h^2)*16*z[,ind[2,i]]^2*z[,ind[1,j]]^2*z[,ind[2,j]]^2)}
			else if (ind[2,i]==ind[1,j]){BB[i,j]=mean((h^2)*16*z[,ind[1,i]]^2*z[,ind[2,j]]^2*z[,ind[1,j]]^2)}
		
		}
	}

	AC=diag(as.vector(C1))



	BC=matrix(0,qind,sp)
	for (i in 1:qind)
	{
		for (j in 1:sp)
		{
			if (j==ind[1,i]){BC[i,j]=mean(h^2*z[,ind[2,i]]^2*z[,ind[1,i]]^2*8)} 
			else if (j==ind[2,i]){BC[i,j]=mean(h^2*z[,ind[1,i]]^2*z[,ind[2,i]]^2*8)} 
		}
	}

	C2=matrix(0,1,p-1)
	for (j in 1:sp)
	{
		C2[j]=mean(4*(h^2)*z[,j]^2)
	}

	CC=diag(as.vector(C2))



	AA2=AA
	for (i in 1:sp)
	{
		for (j in 1:sp)
		{
			AA2[i,j]=16*mean((h^2)*z[,i]^4*z[,j]^4)
		}
	}



	AB2=AB
	for (i in 1:sp)
	{
		for (j in 1:qind)
		{
			AB2[i,j]=32*mean((h^2)*z[,i]^4*z[,ind[1,j]]^2*z[,ind[2,j]]^2)
		}
	}



	BB2=BB
	for (i in 1:qind)
	{
		for (j in 1:qind)
		{
			BB2[i,j]=64*mean((h^2)*(z[,ind[1,i]]^2*z[,ind[2,i]]^2)*(z[,ind[1,j]]^2*z[,ind[2,j]]^2))
		}
	}


	BC2=BC
	for (i in 1:qind)
	{
		for (j in 1:sp)
		{
			BC2[i,j]=16*mean((h^2)*z[,j]^2*(z[,ind[1,i]]^2*z[,ind[2,i]]^2))
		}
	}



	AC2=AC
	for (i in 1:sp)
	{
		for (j in 1:sp)
		{
			AC2[i,j]=8*mean((h^2)*z[,i]^4*z[,j]^2)
		}
	}

	CC2=AC
	for (i in 1:sp)
	{
		for (j in 1:sp)
		{
			CC2[i,j]=4*mean((h^2)*z[,i]^2*z[,j]^2)
		}
	}



	W1=cbind(AA,AB,AC)
	W2=cbind(t(AB),BB,BC)
	W3=cbind(t(AC),t(BC),CC)
	W0=rbind(W1,W2,W3)


	W12=cbind(AA2,AB2,AC2)
	W22=cbind(t(AB2),BB2,BC2)
	W32=cbind(t(AC2),t(BC2),CC2)
	W2=rbind(W12,W22,W32)

	W=W0-W2


	###################
	##calculate d(6)
	##################

	D1=matrix(0,1,sp)
	for (j in 1:sp)
	{
		D1[j]=mean(4*(h^2)*z[,j]^2)
	}
	AD1=diag(as.vector(D1))
	AD=matrix(0,sp,p)
	AD[1:sp,1:sp]=AD1



	BD=matrix(0,qind,p)
	for (i in 1:qind)
	{
		for (j in 1:p)
		{
			if (j==ind[1,i]){BD[i,j]=mean(4*h^2*z[,ind[2,i]]^2)} 
			else if (j==ind[2,i]){BD[i,j]=mean(4*h^2*z[,ind[1,i]]^2)} 
		}
	}


	D2=matrix(0,1,sp)
	for (j in 1:sp)
	{
		D2[j]=mean(2*(h^2))
	}
	CD0=diag(as.vector(D2))
	CD=matrix(0,sp,p)
	CD[1:sp,1:sp]=CD0


	AD2=AD
	for (i in 1:sp)
	{
		for (j in 1:p)
		{
			AD2[i,j]=4*mean((h^2)*z[,i]^4)
		}
	}

	BD2=BD
	for (i in 1:qind)
	{
		for (j in 1:p)
		{
			BD2[i,j]=8*mean((h^2)*(z[,ind[1,i]]^2*z[,ind[2,i]]^2))
		}
	}

	CD2=CD
	for (i in 1:sp)
	{
		for (j in 1:p)
		{
			CD2[i,j]=2*mean((h^2)*z[,i]^2)
		}
	}


	Wnew1=rbind(AD,BD,CD)

	Wnew2=rbind(AD2,BD2,CD2)

	Wnew=Wnew1-Wnew2

	pi2=1+2*beta0

	ev=-1*Wnew%*%pi2
	
	###################
	##calculate d(1)
	##################

	
	lambda4=4*(4+p-2)
	lambda2=2*(2+p-2)

	d_A=matrix(0,sp,1)
	for (j in 1:sp)
	{
		d_A[j]=mean((h^2)*(lambda4*z[,j]^4-12*z[,j]^2))
	}

	d_B=matrix(0,qind,1)
	for (j in 1:qind)
	{
		d_B[j]=2*mean((h^2)*(lambda4*z[,ind[1,j]]^2*z[,ind[2,j]]^2-2*z[,ind[1,j]]^2-2*z[,ind[2,j]]^2))
	}

	d_C=matrix(0,sp,1)
	for (j in 1:sp)
	{
		d_C[j]=mean((h^2)*(lambda2*z[,j]^2-2))
	}

	d=rbind(d_A,d_B,d_C)

	###################
	##calculate d(2)
	##################



	dv_A=matrix(0,n,sp)
	for (i in 1:n)
	{
		for (j in 1:sp)
		{
			if (indh[i]==j){dv_A[i,j]=8*z[i,j]^4*(1-z[i,j]^2)}
			else if (indh[i]==0){dv_A[i,j]=0}
			else {dv_A[i,j]=-8*(z[i,j]^4*z[i,indh[i]]^2)}
		}
	}
	dv_A_mean=matrix(0,1,sp)
	for (j in 1:sp)
	{
		dv_A_mean[j]=mean(dv_A[,j])
	}


	dv_B=matrix(0,n,qind)
	for (i in 1:n)
	{
		for (j in 1:qind)
		{
			if (indh[i]==ind[1,j]){dv_B[i,j]=8*z[i,ind[2,j]]^2*z[i,ind[1,j]]^2*(1-z[i,ind[1,j]]^2)-8*z[i,ind[1,j]]^4*z[i,ind[2,j]]^2}
			else if (indh[i]==ind[2,j]){dv_B[i,j]=8*z[i,ind[2,j]]^2*z[i,ind[1,j]]^2*(1-z[i,ind[2,j]]^2)-8*z[i,ind[1,j]]^2*z[i,ind[2,j]]^4}
			else if (indh[i]==0){dv_B[i,j]=0}
			else {dv_B[i,j]=-16*(z[i,ind[1,j]]^2*z[i,ind[2,j]]^2*z[i,indh[i]]^2)}
		}
	}
	dv_B_mean=matrix(0,1,qind)
	for (j in 1:qind)
	{
		dv_B_mean[j]=mean(dv_B[,j])
	}

	dv_C=matrix(0,n,sp)
	for (i in 1:n)
	{
		for (j in 1:sp)
		{
			if (indh[i]==j){dv_C[i,j]=4*z[i,j]^2*(1-z[i,j]^2)}
			else if (indh[i]==0){dv_C[i,j]=0}
			else {dv_C[i,j]=-4*(z[i,j]^2*z[i,indh[i]]^2)}
		}
	}
	dv_C_mean=matrix(0,1,sp)
	for (j in 1:sp)
	{
		dv_C_mean[j]=mean(dv_C[,j])
	}

	dv=t(cbind(dv_A_mean,dv_B_mean,dv_C_mean))

	###################
	##calculate d total and scoring estimate
	##################

	d=d-dv+ev

	if (incb==1)
	{
		#include bL in the model
		num1=sp+qind+sp
	}
	else
	{
		#omit bL from the model
		num1=sp+qind
	}

	#save W
	W_est=W
	
	#scoring estimator
	quartic_sphere=solve(W[1:num1,1:num1])%*%t(t(d[1:num1]))

	return(list(estimator1=quartic_sphere,W_est=W_est))
}



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
