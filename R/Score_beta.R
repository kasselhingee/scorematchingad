#' @title Score matching estimators for PPI model include beta
#' @description Score matching estimators for the PPI model that estimate \eqn{A_L}, \eqn{b_L} and \eqn{\beta}{beta}.
#' @param prop compositional data (n by p matrix)
#' @param acut \eqn{a_c} for the weighting function \eqn{h}.
#' @param betap the pth element of \eqn{\beta}{beta}, if NULL then this element is estimated.
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
estimatorall1 <- function(prop,acut,betap = NULL)
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


	####calculate W11 ##################
	sp <- p - 1
	ind_qind <- indexcombinations(sp)
	ind <- ind_qind$ind
	qind <- ind_qind$qind
	W11 <- calcW11(p, z, h, ind, qind)

	################### ##calculate W12 ##################
	W12 <- calcW12(p, sp, z, h, ind, qind) #the final column corresponds to the pth element of beta

	################### ##calculate d(1) A ##################
	d1A <- calcd1A(p, sp, z, h, ind, qind)


	################### ##calculate d(2) A ##################
	d2A <- calcd2A_minimah(sp, n, z, ind, qind, indh)

	############################################## ##dirichlet part ###############################################
	Wdir <- calcW22(p, sp, n, z, h, ind, qind)
	h4m <- h2onz2_mean(p, n, z, h, indh)

	##### Calculate d(1) B #####
	d1=t(((p-2)*mean(h^2)+h4m))  #this is d(1)_B as defined in Section A.5

  ##### Calculate d(2) B #####
	# is sensitive to choise of weight function
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

  ################### Calculate d(6) when beta0[p] is fixed ######
	# this is -W12 * pi for the fixed beta, For this function beta[p] may be fixed
	if (!is.null(betap)){
	  pi2 <- rep(0, p)
	  pi2[p] <- 1 + 2*betap
	  d6=-1*rbind(W12, Wdir) %*%pi2
	}

	################### ##calculate d total and scoring estimate ##################

	dA=d1A+d2A

	d=rbind(dA,ddir)

	W=rbind(cbind(W11,W12),cbind(t(W12),Wdir))

	if (is.null(betap))
	{
		#estimate beta[p] in the model
		num1=sp+qind+sp+p
	}
	else
	{
		#fix beta[p] in model
		num1=sp+qind+sp+sp
		# add offset
		d <- d + d6
	}

	#save W
	W_est=W

	#scoring estimator
	quartic_sphere=solve(W[1:num1,1:num1])%*%t(t(d[1:num1]))

	tot=sp+qind+sp

	# from estimates of 1 + 2beta, to estimates of beta
	quartic_sphere[sum(tot,1):num1]=(quartic_sphere[sum(tot,1):num1]-1)/2

	return(list(estimator1=quartic_sphere,W_est=W_est))
}



