#' @title Score matching estimate of the PPI model using a minima-based Hyvarinen weight function
#' @description Estimates \eqn{A_L} and \eqn{b_L} of the PPI model using score matching and a minimum-like Hyvarinen weight function. \eqn{\beta_0}{beta0} is fixed.
#' @param prop compositional data (each row is a sample, each column corresponds to a component)
#' @param  acut \eqn{a_c} for the weighting function \eqn{h}.
#' @param incb if `incb=1` then \eqn{b_L} is estimated, otherwise \eqn{b_L} is fixed at zero
#' @param beta0 The (fixed) beta0 of the model.
#' @details The PPI model is given in equation 3 of (Scealy and Wood, 2021). The matrices \eqn{A_L} and \eqn{b_L} must be estimated.
#' This function implements the score matching estimator,
#' \deqn{\hat{W}^{-1}\hat{d},}{W^{-1}d,}
#' using a minima-based Hyvarinen weight function
#' \deqn{\tilde{h}(z)^2 = \min(z_1^2, z_2^2, ..., z_p^2, a_c^2).}{h(z)^2 = min(z1^2, z2^2, ..., zp^2, a_c^2),}
#' where \eqn{z} is a point in the positive orthant of the p-dimensional unit sphere
#' and \eqn{z_j}{zj} is the jth component of z.
#' For details of the score matching estimator see equations 16 - 19 in (Scealy and Wood, 2021).
#' If \eqn{a_c} is greater than or equal to 1 then this Hyvarinen weight function corresponds to (Scealy and Wood, 2021; eqn 11), if it is less than 1 then it corresponds to (Scealy and Wood, 2021; eqn 12).
#' For more on the Hyvarinen weight (see equation 7 and Section 3.2 of (Scealy and Wood, 2021)).
#' @return A vector of the estimates for individual entries in the matrices \eqn{A} and \eqn{b}, and the estimated \eqn{\hat{W}}{W}. The former first contains the diagonal of \eqn{A} (except the last entry that is always zero for identifiability in the PPI model), then the upper triangle of \eqn{A} without the last column (again for identifiability) and finally the elements of \eqn{b} (except the last element, which is always 0 due to identifiability also) if `incb=1`.


#' @export
estimator1 <- function(prop,acut,incb, beta0)
{
  n <- nrow(prop) #number of samples
  p <- ncol(prop) #number of dimensions, although what happens when the beta need to be estimated?
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


