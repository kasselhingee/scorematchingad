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


	################### ##calculate W ##################
	sp <- p - 1
	ind_qind <- indexcombinations(sp)
	ind <- ind_qind$ind
	qind <- ind_qind$qind
	W <- calcW11(p, z, h, ind, qind)

	################### ##calculate d(6) ##################
	ev <- calcd6(p, sp, z, h, ind, qind, beta0)

	################### ##calculate d(1) ##################
	d <- calcd1A(p, sp, z, h, ind, qind)

	################### ##calculate d(2) ##################
	dv <- calcd2A_minimah(sp, n, z, ind, qind, indh)

	################### ##calculate d total and scoring estimate ##################

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


