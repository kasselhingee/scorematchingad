#' @title Score matching estimators for PPI model include beta
#' @description Score matching estimators for the PPI model that estimate \eqn{A_L}, \eqn{b_L} and \eqn{\beta}{beta}.
#' @param prop compositional data (n by p matrix)
#' @param acut \eqn{a_c} for the weighting function \eqn{h}.
#' @param incb if `incb=1` then the pth element of \eqn{\beta}{beta} is estimated, otherwise it is fixed at -0.5.
#' @details The PPI model is given in equation 3 of (Scealy and Wood, 2021). The matrix \eqn{A_L} and vectors \eqn{b_L} and \eqn{\beta}{beta} must be estimated.
#' This function implements the score matching estimator,
#' \deqn{\hat{W}^{-1}\hat{d},}{W^{-1}d,}
#' using a minima-based Hyvarinen weight function
#' \deqn{\tilde{h}(z)^2 = \min(z_1^2, z_2^2, ..., z_p^2, a_c^2).}{h(z)^2 = min(z1^2, z2^2, ..., zp^2, a_c^2),}
#' where \eqn{z} is a point in the positive orthant of the p-dimensional unit sphere
#' and \eqn{z_j}{zj} is the jth component of z.
#' For details of the score matching estimator see equations 16 - 19 in (Scealy and Wood, 2021).
#' If \eqn{a_c} is greater than or equal to 1 then this Hyvarinen weight function corresponds to (Scealy and Wood, 2021; eqn 11), if it is less than 1 then it corresponds to (Scealy and Wood, 2021; eqn 12).
#' For more on the Hyvarinen weight (see equation 7 and Section 3.2 of (Scealy and Wood, 2021)).

#' @return A vector of the estimates for individual entries of \eqn{A_L}, \eqn{b_L}, and \eqn{\beta}{beta}, and the estimated \eqn{\hat{W}}{W}. The former first contains the diagonal of \eqn{A_L}, then the upper triangle of \eqn{A_L}, then the elements of \eqn{b_L}, and then finally the estimates of \eqn{\beta}{beta}.

#' @export
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


	####calculate W ##################
	sp <- p - 1
	ind_qind <- indexcombinations(sp)
	ind <- ind_qind$ind
	qind <- ind_qind$qind
	W <- calcW11(p, z, h, ind, qind)

	Whybrid=W

	################### ##calculate d(6) ##################
	Wboth <- calcW12(p, sp, z, h, ind, qind) #the final column corresponds to the pth element of beta

	#pi2=1+2*beta0

	#ev=-1*Wnew%*%pi2

	################### ##calculate d(1) ##################
	d <- calcd1A(p, sp, z, h, ind, qind)

	d1hybrid=d

	################### ##calculate d(2) ##################
	dv <- calcd2A_minimah(sp, n, z, ind, qind, indh)
	d2hybrid=dv

	############################################## ##dirichlet part ###############################################

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
	ind=utils::combn(x, 2, FUN = NULL, simplify = TRUE)
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



	################### ##calculate d total and scoring estimate ##################

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



