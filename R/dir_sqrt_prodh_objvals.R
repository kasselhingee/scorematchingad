#' @title Value of, up to a constant, components of score matching objective for the Dirichlet distribution
#' @description Computes the matrices of the score matching objective in eqn 15 for the Dirichlet distribution, which is the PPI with \eqn{A_L=0} and \eqn{b_L=0}. The parameters to be estimates are the \eqn{\beta}{beta}.
#' @param dirfit Compositional data (n by p matrix)
#' @param acut \eqn{a_c} for \eqn{h} function.
#' @return A list of matrices
#' \deqn{\tilde{h}(z)^2 = \min(\prod_{j=1}^{p} z_j^2, a_c^2).}{h(z)^2 = min(z1^2 * z2^2 * ... * zp^2, a_c^2).}
#' @export
score2_mats_dir <- function(dirfit,acut)
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



	sp=p-1

        indqind <- indexcombinations(sp)
	ind=indqind$ind
	qind=indqind$qind

	h4m <- h2onz2_mean(p, n, z, h, indh, hstyle = "product")
  W <- calcW22(p, h, h4m)

	d1=t(((p-2)*mean(h^2)+h4m))

	#ind2 indicates whether the acut constraint was hit
	ind2=matrix(1,n,1)
	for (j in 1:n)
	{
		if (h[j] == acut) {ind2[j]=0}
	}

  # the following is h4m as above, except that h4s is zero
	# whenever the acut constraint is hit
	# this is the first half of the equation for d(2)B at the end of
	# page 22 in the appendix
	#homit matrix. For a given component, k, homit is the multiple
	# of the z entries in the other components. I.e. prod(z_j)/z_k when z_k is not zero.
	homit=matrix(1,n,p)
	for (k in 1:p)
	{
		for (j in 1:p)
		{
			if (k != j){homit[,k]=homit[,k]*z[,j]}
		}

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

	out <- list(
	  W = W,
	  d1 = d1,
	  d2 = -d2,
	  h = h
	)

	return(out)
}

#' @describeIn score2_mats_dir The score matching estimator using the product-based Hyvarinen weight
#' @export
score2_dir <- function(dirfit, acut, beta){
  mats <- score2_mats_dir(dirfit, acut)
  pi <- matrix(1 + 2*beta, ncol = 1)
  out <- c(
    hgradpgrad = 0.5 * t(pi) %*% mats$W %*% pi,
    hlap = - t(pi) %*% mats$d1,
    gradhPgrad = -t(pi) %*% mats$d2
  )
  return(out)
}
