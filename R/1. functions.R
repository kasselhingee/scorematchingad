######################################################################
##rhybrid generates random observations from the PPI model
##n: sample size
##p: dimension
##beta0: beta shape parameter vector
##ALs: A_L parameter matrix
##bL: b_L parameter vector
##maxden: is the constant log(C) in Appendix A.1.3.
######################################################################

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


#############################################################
##rtGaussian generates random observations from the truncated
##Gaussian model (beta=0)
##n: sample size
##p: dimension
##muL: untruncated Gaussian mean (see Appendix A.1.1)
##SigA: untruncated Gaussian covariance matrix (see Appendix A.1.1)
#############################################################


library(MASS)

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





#######################################################
##estimator1SE calculates Score1ac standard errors in the article
##for the model with A_L and b_L estimated and beta fixed.
##prop: compositional data (n by p matrix)
##acut: a_c for h function
##estimate1: the value of the Score1ac estimate
##W_est: Estimated W matrix
##incb: if incb=1 then b_L is estimated otherwise b_L is fixed at zero (omitted)
##beta0: passed beta0
#######################################################




estimator1SE <- function(prop,acut,estimate1,W_est,incb, beta0)
{
  n<-nrow(prop)
  p<-ncol(prop)


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


	sp=p-1
	x=c(1:sp)
	ind=combn(x, 2, FUN = NULL, simplify = TRUE)
	qind=length(ind[1,])

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



	AD_big=array(0, dim=c(sp, p, n))
	for (i in 1:n)
	{
		D1=matrix(0,1,sp)
		for (j in 1:sp)
		{
			D1[j]=4*(h[i]^2)*z[i,j]^2
		}
		AD1=diag(as.vector(D1))
		AD=matrix(0,sp,p)
		AD[1:sp,1:sp]=AD1
		AD_big[1:sp,1:p,i]=AD
	}


	BD_big=array(0, dim=c(qind, p, n))
	for (k in 1:n)
	{
		BD=matrix(0,qind,p)
		for (i in 1:qind)
		{
			for (j in 1:p)
			{
				if (j==ind[1,i]){BD[i,j]=mean(4*h[k]^2*z[k,ind[2,i]]^2)}
				else if (j==ind[2,i]){BD[i,j]=mean(4*h[k]^2*z[k,ind[1,i]]^2)}
			}
		}
		BD_big[1:qind,1:p,k]=BD
	}


	CD_big=array(0, dim=c(sp, p, n))
	for (k in 1:n)
	{
		D2=matrix(0,1,sp)
		for (j in 1:sp)
		{
			D2[j]=mean(2*(h[k]^2))
		}
		CD0=diag(as.vector(D2))
		CD=matrix(0,sp,p)
		CD[1:sp,1:sp]=CD0
		CD_big[1:sp,1:p,k]=CD
	}


	AD2_big=array(0, dim=c(sp, p, n))
	AD2=AD
	for (k in 1:n)
	{
		for (i in 1:sp)
		{
			for (j in 1:p)
			{
				AD2[i,j]=4*mean((h[k]^2)*z[k,i]^4)
			}
		}
		AD2_big[1:sp,1:p,k]=AD2
	}



	BD2_big=array(0, dim=c(qind, p, n))
	BD2=BD
	for (k in 1:n)
	{
		for (i in 1:qind)
		{
			for (j in 1:p)
			{
				BD2[i,j]=8*mean((h[k]^2)*(z[k,ind[1,i]]^2*z[k,ind[2,i]]^2))
			}
		}
		BD2_big[1:qind,1:p,k]=BD2
	}


	CD2_big=array(0, dim=c(sp, p, n))
	CD2=CD
	for (k in 1:n)
	{
		for (i in 1:sp)
		{
			for (j in 1:p)
			{
				CD2[i,j]=2*mean((h[k]^2)*z[k,i]^2)
			}
		}
		CD2_big[1:sp,1:p,k]=CD2
	}

	pi2=1+2*beta0
	evm=matrix(0,sum(sp,sp,qind),n)
	for (k in 1:n)
	{

		Wnew1=rbind(AD_big[1:sp,1:p,k],BD_big[1:qind,1:p,k],CD_big[1:sp,1:p,k])
		Wnew2=rbind(AD2_big[1:sp,1:p,k],BD2_big[1:qind,1:p,k],CD2_big[1:sp,1:p,k])
		Wnew=Wnew1-Wnew2


		evm[,k]=-1*Wnew%*%pi2
	}



	lambda4=4*(4+p-2)
	lambda2=2*(2+p-2)

	d_A=matrix(0,sp,n)
	for (j in 1:sp)
	{
		d_A[j,]=((h^2)*(lambda4*z[,j]^4-12*z[,j]^2))
	}

	d_B=matrix(0,qind,n)
	for (j in 1:qind)
	{
		d_B[j,]=2*((h^2)*(lambda4*z[,ind[1,j]]^2*z[,ind[2,j]]^2-2*z[,ind[1,j]]^2-2*z[,ind[2,j]]^2))
	}

	d_C=matrix(0,sp,n)
	for (j in 1:sp)
	{
		d_C[j,]=((h^2)*(lambda2*z[,j]^2-2))
	}

	dm=rbind(d_A,d_B,d_C)

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


	dvm=t(cbind(dv_A,dv_B,dv_C))

	dm=dm-dvm+evm


	AA_big=array(0, dim=c(sp, sp, n))
	for (k in 1:n)
	{
		A1=matrix(0,1,p-1)
		for (j in 1:sp)
		{
			A1[j]=(16*(h[k]^2)*z[k,j]^6)
		}
		AA=diag(as.vector(A1))
		AA_big[1:sp,1:sp,k]=AA
	}



	AB_big=array(0, dim=c(sp, qind, n))
	for (k in 1:n)
	{
		AB=matrix(0,sp,qind)
		for (i in 1:sp)
		{
			for (j in 1:qind)
			{
				if (i==ind[1,j]){AB[i,j]=(16*z[k,ind[1,j]]^4*z[k,ind[2,j]]^2*h[k]^2)}
				if (i==ind[2,j]){AB[i,j]=(16*z[k,ind[1,j]]^2*z[k,ind[2,j]]^4*h[k]^2)}
			}
		}
		AB_big[1:sp,1:qind,k]=AB
	}


	BB_big=array(0, dim=c(qind, qind, n))
	for (k in 1:n)
	{
		BB=matrix(0,qind,qind)
		for (i in 1:qind)
		{
			for (j in 1:qind)
			{
				if (ind[1,i]==ind[1,j] & ind[2,i]==ind[2,j]){BB[i,j]=(16*z[k,ind[1,j]]^4*z[k,ind[2,j]]^2*h[k]^2)+(16*z[k,ind[1,j]]^2*z[k,ind[2,j]]^4*h[k]^2)}
				else if (ind[1,i]==ind[1,j]){BB[i,j]=((h[k]^2)*16*z[k,ind[2,j]]^2*z[k,ind[2,i]]^2*z[k,ind[1,j]]^2)}
				else if (ind[2,i]==ind[2,j]){BB[i,j]=((h[k]^2)*16*z[k,ind[1,j]]^2*z[k,ind[1,i]]^2*z[k,ind[2,j]]^2)}
				else if (ind[1,i]==ind[2,j]){BB[i,j]=((h[k]^2)*16*z[k,ind[2,i]]^2*z[k,ind[1,j]]^2*z[k,ind[2,j]]^2)}
				else if (ind[2,i]==ind[1,j]){BB[i,j]=((h[k]^2)*16*z[k,ind[1,i]]^2*z[k,ind[2,j]]^2*z[k,ind[1,j]]^2)}

			}
		}
		BB_big[1:qind,1:qind,k]=BB
	}



	AC_big=array(0, dim=c(sp, sp, n))
	for (k in 1:n)
	{

		C1=matrix(0,1,p-1)
		for (j in 1:sp)
		{
			C1[j]=(8*(h[k]^2)*z[k,j]^4)
		}
		AC=diag(as.vector(C1))
		AC_big[1:sp,1:sp,k]=AC
	}




	BC_big=array(0, dim=c(qind, sp, n))
	for (k in 1:n)
	{

		BC=matrix(0,qind,sp)
		for (i in 1:qind)
		{
			for (j in 1:sp)
			{
				if (j==ind[1,i]){BC[i,j]=(h[k]^2*z[k,ind[2,i]]^2*z[k,ind[1,i]]^2*8)}
				else if (j==ind[2,i]){BC[i,j]=(h[k]^2*z[k,ind[1,i]]^2*z[k,ind[2,i]]^2*8)}
			}
		}
		BC_big[1:qind,1:sp,k]=BC
	}


	CC_big=array(0, dim=c(sp, sp, n))
	for (k in 1:n)
	{

		C2=matrix(0,1,p-1)
		for (j in 1:sp)
		{
			C2[j]=(4*(h[k]^2)*z[k,j]^2)
		}
		CC=diag(as.vector(C2))
		CC_big[1:sp,1:sp,k]=CC
	}


	AA2_big=array(0, dim=c(sp, sp, n))
	for (k in 1:n)
	{

		AA2=AA
		for (i in 1:sp)
		{
			for (j in 1:sp)
			{
				AA2[i,j]=16*((h[k]^2)*z[k,i]^4*z[k,j]^4)
			}
		}
		AA2_big[1:sp,1:sp,k]=AA2
	}



	AB2_big=array(0, dim=c(sp, qind, n))
	for (k in 1:n)
	{

		AB2=AB
		for (i in 1:sp)
		{
			for (j in 1:qind)
			{
				AB2[i,j]=32*((h[k]^2)*z[k,i]^4*z[k,ind[1,j]]^2*z[k,ind[2,j]]^2)
			}
		}
		AB2_big[1:sp,1:qind,k]=AB2
	}


	BB2_big=array(0, dim=c(qind, qind, n))
	for (k in 1:n)
	{
		BB2=BB
		for (i in 1:qind)
		{
			for (j in 1:qind)
			{
				BB2[i,j]=64*((h[k]^2)*(z[k,ind[1,i]]^2*z[k,ind[2,i]]^2)*(z[k,ind[1,j]]^2*z[k,ind[2,j]]^2))
			}
		}
		BB2_big[1:qind,1:qind,k]=BB2
	}



	BC2_big=array(0, dim=c(qind, sp, n))
	for (k in 1:n)
	{

		BC2=BC
		for (i in 1:qind)
		{
			for (j in 1:sp)
			{
				BC2[i,j]=16*((h[k]^2)*z[k,j]^2*(z[k,ind[1,i]]^2*z[k,ind[2,i]]^2))
			}
		}
		BC2_big[1:qind,1:sp,k]=BC2
	}


	AC2_big=array(0, dim=c(sp, sp, n))
	for (k in 1:n)
	{
		AC2=AC
		for (i in 1:sp)
		{
			for (j in 1:sp)
			{
				AC2[i,j]=8*((h[k]^2)*z[k,i]^4*z[k,j]^2)
			}
		}
		AC2_big[1:sp,1:sp,k]=AC2
	}



	CC2_big=array(0, dim=c(sp, sp, n))
	for (k in 1:n)
	{
		CC2=AC
		for (i in 1:sp)
		{
			for (j in 1:sp)
			{
				CC2[i,j]=4*((h[k]^2)*z[k,i]^2*z[k,j]^2)
			}
		}
		CC2_big[1:sp,1:sp,k]=CC2
	}


	diff=matrix(0,sum(sp,sp,qind),n)
	for (k in 1:n)
	{

		AA[,]=AA_big[,,k]
		AB[,]=AB_big[,,k]
		AC[,]=AC_big[,,k]

		AB[,]=AB_big[,,k]
		BB[,]=BB_big[,,k]
		BC[,]=BC_big[,,k]

		AC[,]=AC_big[,,k]
		BC[,]=BC_big[,,k]
		CC[,]=CC_big[,,k]


		W1=cbind(AA,AB,AC)
		W2=cbind(t(AB),BB,BC)
		W3=cbind(t(AC),t(BC),CC)
		W0=rbind(W1,W2,W3)

		AA2[,]=AA2_big[,,k]
		AB2[,]=AB2_big[,,k]
		AC2[,]=AC2_big[,,k]

		AB2[,]=AB2_big[,,k]
		BB2[,]=BB2_big[,,k]
		BC2[,]=BC2_big[,,k]

		AC2[,]=AC2_big[,,k]
		BC2[,]=BC2_big[,,k]
		CC2[,]=CC2_big[,,k]


		W12=cbind(AA2,AB2,AC2)
		W22=cbind(t(AB2),BB2,BC2)
		W32=cbind(t(AC2),t(BC2),CC2)
		W2=rbind(W12,W22,W32)

		W=W0-W2

		diff[1:num1,k]=t(t(dm[1:num1,k]))-W[1:num1,1:num1]%*%estimate1
	}

	Sig0=(diff[1:num1,]%*%t(diff[1:num1,]))/n
	Gam0=W_est[1:num1,1:num1]
	var0=solve(Gam0)%*%Sig0%*%solve(Gam0)/n


	std=sqrt(diag(var0))



	return(std)
}


#######################################################
##multestimator calculates the ScoreMult estimator in the article (p=3 only)
##and assumes beta=beta0 is the same in each category.
##Note: beta0 is fixed and not estimated.
##Currently only A_L is estimated (b_L is fixed at zero, but this can
##be modified by including components 4 and 5 in the estimator
##pp in the last part of the code).
##x: count data (n by 3 matrix)
##ni: see 5. model1.R for an example of use
#######################################################




multestimator <- function(x, ni, beta0)
{

	X1=x
	X2=x*(x-1)
	X3=x*(x-1)*(x-2)
	X4=x*(x-1)*(x-2)*(x-3)
	X5=x*(x-1)*(x-2)*(x-3)*(x-4)
	X6=x*(x-1)*(x-2)*(x-3)*(x-4)*(x-5)

	fac1=ni
	fac2=ni*(ni-1)
	fac3=ni*(ni-1)*(ni-2)
	fac4=ni*(ni-1)*(ni-2)*(ni-3)
	fac5=ni*(ni-1)*(ni-2)*(ni-3)*(ni-4)
	fac6=ni*(ni-1)*(ni-2)*(ni-3)*(ni-4)*(ni-5)
	fac7=ni*(ni-1)*(ni-2)*(ni-3)*(ni-4)*(ni-5)*(ni-6)

	p3p1=mean(X3[,1]*X1[,2]/fac4)
	p4p1=mean(X4[,1]*X1[,2]/fac5)
	p3p2=mean(X3[,1]*X2[,2]/fac5)
	p2p1=mean(X2[,1]*X1[,2]/fac3)
	p2p2=mean(X2[,1]*X2[,2]/fac4)
	p1p3=mean(X1[,1]*X3[,2]/fac4)
	p1p4=mean(X1[,1]*X4[,2]/fac5)
	p2p3=mean(X2[,1]*X3[,2]/fac5)
	p1p2=mean(X1[,1]*X2[,2]/fac3)
	p2p2=mean(X2[,1]*X2[,2]/fac4)
	p3p3=mean(X3[,1]*X3[,2]/fac6)
	p1p1=mean(X1[,1]*X1[,2]/fac2)
	p5p1=mean(X5[,1]*X1[,2]/fac6)
	p4p2=mean(X4[,1]*X2[,2]/fac6)
	p6p1=mean(X6[,1]*X1[,2]/fac7)
	p5p2=mean(X5[,1]*X2[,2]/fac7)
	p3p4=mean(X3[,1]*X4[,2]/fac7)
	p4p3=mean(X4[,1]*X3[,2]/fac7)
	p2p4=mean(X2[,1]*X4[,2]/fac6)
	p1p5=mean(X1[,1]*X5[,2]/fac6)
	p2p5=mean(X2[,1]*X5[,2]/fac7)
	p1p6=mean(X1[,1]*X6[,2]/fac7)



	d=matrix(0,5,1)
	d[1]=32*p3p1-20*p4p1-20*p3p2-12*p2p1+12*p2p2
	d[2]=32*p1p3-20*p1p4-20*p2p3-12*p1p2+12*p2p2
	d[3]=48*p2p2-40*p3p2-40*p2p3-4*p1p2+4*p1p3-4*p2p1+4*p3p1
	d[4]=8*p2p1-6*p3p1-2*p1p1+2*p1p2-6*p2p2
	d[5]=8*p1p2-6*p2p2-6*p1p3-2*p1p1+2*p2p1



	dv=matrix(0,5,1)
	dv[1]=2*(4*p2p1-16*p3p1-4*p2p2+12*p4p1+12*p3p2)
	dv[2]=2*(4*p1p2-4*p2p2-16*p1p3+12*p2p3+12*p1p4)
	dv[3]=2*(4*p1p2+4*p2p1-32*p2p2-4*p1p3-4*p3p1+24*p3p2+24*p2p3)
	dv[4]=2*(2*p1p1-8*p2p1-2*p1p2+6*p3p1+6*p2p2)
	dv[5]=2*(2*p1p1-8*p1p2-2*p2p1+6*p1p3+6*p2p2)


	ev=matrix(0,5,1)
	ev[1]=12*(p3p1-p4p1-p3p2)-4*(p2p1-p3p1-p2p2)
	ev[2]=12*(p1p3-p2p3-p1p4)-4*(p1p2-p2p2-p1p3)
	ev[3]=24*(p2p2-p3p2-p2p3)-4*(p2p1-p3p1-p2p2)-4*(p1p2-p2p2-p1p3)
	ev[4]=6*(p2p1-p3p1-p2p2)-2*(p1p1-p2p1-p1p2)
	ev[5]=6*(p1p2-p2p2-p1p3)-2*(p1p1-p2p1-p1p2)



	d=d-dv+ev*(1+2*beta0[1])



	W=matrix(0,5,5)
	W[1,1]=16*p4p1-32*p5p1-16*p4p2+16*p6p1+16*p5p2
	W[1,2]=-16*p3p3+16*p3p4+16*p4p3
	W[1,3]=16*p3p2-48*p4p2-16*p3p3+32*p5p2+32*p4p3
	W[1,4]=8*p3p1-16*p4p1-8*p3p2+8*p5p1+8*p4p2
	W[1,5]=-8*p3p2+8*p4p2+8*p3p3
	W[2,1]=W[1,2]
	W[3,1]=W[1,3]
	W[4,1]=W[1,4]
	W[5,1]=W[1,5]
	W[2,2]=16*p1p4-16*p2p4-32*p1p5+16*p2p5+16*p1p6
	W[2,3]=16*p2p3-16*p3p3-48*p2p4+32*p3p4+32*p2p5
	W[2,4]=-8*p2p3+8*p3p3+8*p2p4
	W[2,5]=8*p1p3-8*p2p3-16*p1p4+8*p2p4+8*p1p5
	W[3,2]=W[2,3]
	W[4,2]=W[2,4]
	W[5,2]=W[2,5]
	W[3,3]=16*p2p3-96*p3p3-16*p2p4-16*p4p2+16*p3p2+64*p4p3+64*p3p4
	W[3,4]=8*p2p2-24*p3p2-8*p2p3+16*p4p2+16*p3p3
	W[3,5]=8*p2p2-8*p3p2-24*p2p3+16*p3p3+16*p2p4
	W[4,3]=W[3,4]
	W[5,3]=W[3,5]
	W[4,4]=4*p2p1-8*p3p1-4*p2p2+4*p4p1+4*p3p2
	W[4,5]=-4*p2p2+4*p3p2+4*p2p3
	W[5,4]=W[4,5]
	W[5,5]=4*p1p2-4*p2p2-8*p1p3+4*p2p3+4*p1p4

	pp=solve(W[1:3,1:3])%*%t(t(d[1:3]))

	mult=pp

	return(mult)

}




#######################################################
##estimator2_dir calculates estimators Score2ac and Score2 in the article
##for the Dirichlet distribution
##dirfit: compositional data (n by p matrix)
##acut: a_c for h function (if acut > 1 then this is Score2)
#######################################################




estimator2_dir <- function(dirfit,acut)
{
  n=nrow(dirfit)
  p=ncol(dirfit)
	z=sqrt(dirfit)

	h=matrix(1,n,1)
	for (j in 1:p)
	{
		h=h*z[,j]
	}
	for (j in 1:n)
	{
		if (h[j] > acut) {h[j]=acut}
	}

	homit=matrix(1,n,p)
	for (k in 1:p)
	{
		for (j in 1:p)
		{
			if (k != j){homit[,k]=homit[,k]*z[,j]}
		}

	}


	sp=p-1

	x=c(1:sp)
	ind=combn(x, 2, FUN = NULL, simplify = TRUE)
	qind=length(ind[1,])

	h4s=matrix(0,n,p)
	h4m=matrix(0,1,p)
	for (j in 1:p)
	{
		for (i in 1:n)
		{
			if (h[i] > 0){h4s[i,j]=(h[i]^2)/z[i,j]^2}
			else {h4s[i,j]=homit[i,j]^2}

		}
		h4m[j]=mean(h4s[,j])
	}


	DD=matrix(0,p,p)
	for (i in 1:p)
	{
		for (j in 1:p)
		{
			if (i==j){DD[i,j]=h4m[i]}
		}
	}

	DD2=matrix(mean(h^2),p,p)


	W=DD-DD2

	d1=t(((p-2)*mean(h^2)+h4m))

	ind2=matrix(1,n,1)
	for (j in 1:n)
	{
		if (h[j] == acut) {ind2[j]=0}
	}


	h4s=matrix(0,n,p)
	h4m=matrix(0,1,p)
	for (j in 1:p)
	{
		for (i in 1:n)
		{
			if (h[i] > 0){h4s[i,j]=ind2[i]*(h[i]^2)/z[i,j]^2}
			else {h4s[i,j]=ind2[i]*homit[i,j]^2}

		}
		h4m[j]=mean(h4s[,j])
	}



	d2=t(2*h4m-2*p*mean(ind2*h^2))

	d=d1-d2

	pp=solve(W)%*%d


	estimate2=(pp-1)/2

	return(estimate2)

}


#######################################################
##estimator1_dir calculates estimators Score1ac and Score1 in the article
##for the Dirichlet distribution
##dirfit: compositional data (n by p matrix)
##acut: a_c for h function (if acut > 1 then this is Score1)
#######################################################



estimator1_dir <- function(dirfit,acut)
{
  n=nrow(dirfit)
  p=ncol(dirfit)

	z=sqrt(dirfit)

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

	sp=p-1

	x=c(1:sp)
	ind=combn(x, 2, FUN = NULL, simplify = TRUE)
	qind=length(ind[1,])

	h4s=matrix(0,n,p)
	h4m=matrix(0,1,p)
	for (j in 1:p)
	{
		for (i in 1:n)
		{
			if (h[i] > 0){h4s[i,j]=(h[i]^2)/z[i,j]^2}
			else if (indh[i]==j){h4s[i,j]=1}

		}
		h4m[j]=mean(h4s[,j])
	}


	DD=matrix(0,p,p)
	for (i in 1:p)
	{
		for (j in 1:p)
		{
			if (i==j){DD[i,j]=h4m[i]}
		}
	}

	DD2=matrix(mean(h^2),p,p)


	W=DD-DD2

	d1=t(((p-2)*mean(h^2)+h4m))


	d2=matrix(0,n,p)

	for (i in 1:n)
	{
		for (j in 1:p)
		{
			if (indh[i]==j){d2[i,j]=2*(1-z[i,j]^2)}
			else if (indh[i]==0){d2[i,j]=0}
			else {d2[i,j]=-2*(z[i,indh[i]]^2)}
		}
	}


	d3=matrix(0,1,p)
	for (j in 1:p)
	{
		d3[j]=mean(d2[,j])
	}


	d=d1-t(d3)

	pp=solve(W)%*%d
	estimate1=(pp-1)/2

	return(estimate1)
}



####################################################################
##method of moments estimator for Dirichlet distribution
####################################################################

dirichmom <- function(X) {
# Method of Moments estimates of the Dirichlet Distribution
temp <- dim(X); n <- temp[1]; m <- temp[2]
# X <- cbind(X,matrix(1-apply(X,1,sum)))
mom <- apply(X,2,mean)*(mean(X[,1])-mean(X[,1]^2))/(mean(X[,1]^2) - ((mean(X[,1]))^2))
return(matrix(mom))
}





#######################################################
##estimatorall1 calculates estimators Score1 and Score1ac in the article
##for the full model (including A_L, b_L and beta all estimated).
##prop: compositional data (n by p matrix)
##acut: a_c for h function (if acut > 1 then this is Score1)
##incb: if incb=1 then beta[p] is estimated otherwise beta[p] is fixed at -0.5
#######################################################



estimatorall1 <- function(prop,acut,incb)
{
  n = nrow(prop)
  p = ncol(prop)

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

	Whybrid=W


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

	Wboth=Wnew

	#pi2=1+2*beta0

	#ev=-1*Wnew%*%pi2

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

	d1hybrid=d

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

	d2hybrid=dv

	##############################################
	##dirichlet part
	###############################################

	z=sqrt(prop)

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

	sp=p-1

	x=c(1:sp)
	ind=combn(x, 2, FUN = NULL, simplify = TRUE)
	qind=length(ind[1,])

	h4s=matrix(0,n,p)
	h4m=matrix(0,1,p)
	for (j in 1:p)
	{
		for (i in 1:n)
		{
			if (h[i] > 0){h4s[i,j]=(h[i]^2)/z[i,j]^2}
			else if (indh[i]==j){h4s[i,j]=1}

		}
		h4m[j]=mean(h4s[,j])
	}



	DD=matrix(0,p,p)
	for (i in 1:p)
	{
		for (j in 1:p)
		{
			if (i==j){DD[i,j]=h4m[i]}
		}
	}

	DD2=matrix(mean(h^2),p,p)


	W=DD-DD2

	Wdir=W

	d1=t(((p-2)*mean(h^2)+h4m))


	d2=matrix(0,n,p)

	for (i in 1:n)
	{
		for (j in 1:p)
		{
			if (indh[i]==j){d2[i,j]=2*(1-z[i,j]^2)}
			else if (indh[i]==0){d2[i,j]=0}
			else {d2[i,j]=-2*(z[i,indh[i]]^2)}
		}
	}


	d3=matrix(0,1,p)
	for (j in 1:p)
	{
		d3[j]=mean(d2[,j])
	}

	d=d1-t(d3)

	ddir=d



	###################
	##calculate d total and scoring estimate
	##################

	dhybrid=d1hybrid-d2hybrid

	d=rbind(dhybrid,ddir)

	W=rbind(cbind(Whybrid,Wboth),cbind(t(Wboth),Wdir))

	if (incb==1)
	{
		#include beta[p] in the model
		num1=sp+qind+sp+p
	}
	else
	{
		#fix beta[p] in model
		num1=sp+qind+sp+sp
	}

	#save W
	W_est=W

	#scoring estimator
	quartic_sphere=solve(W[1:num1,1:num1])%*%t(t(d[1:num1]))

	tot=sp+qind+sp

	quartic_sphere[sum(tot,1):num1]=(quartic_sphere[sum(tot,1):num1]-1)/2

	return(list(estimator1=quartic_sphere,W_est=W_est))
}


