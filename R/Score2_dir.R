#' @title Score matching estimates for the Dirichlet distribution
#' @description Score matching estimates for the Dirichlet distribution, which is the PPI with \eqn{A_L=0} and \eqn{b_L=0}. The parameters to be estimates are the \eqn{\beta}{beta}.
#' @param dirfit Compositional data (n by p matrix)
#' @param acut \eqn{a_c} for \eqn{h} function.
#' @return A vector of estimates.
#' @describeIn estimator_dir The score matching estimator using the product-based Hyvarinen weight
#' \deqn{\tilde{h}(z)^2 = \min(\prod_{j=1}^{p} z_j^2, a_c^2).}{h(z)^2 = min(z1^2 * z2^2 * ... * zp^2, a_c^2).}
#' @export
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

#' @describeIn estimator_dir The score matching estimator using the minima-based Hyvarinen weight function 
#' \deqn{\tilde{h}(z)^2 = \min(z_1^2, z_2^2, ..., z_p^2, a_c^2).}{h(z)^2 = min(z1^2, z2^2, ..., zp^2, a_c^2).}
#' @export
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

