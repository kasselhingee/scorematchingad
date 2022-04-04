
calcW22 <- function(p, sp, n, z, h, ind, qind, indh){
  h4m <- h2onz2_mean(p, n, z, h, indh)

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
  return(W)
}

# compute the first half of of the sample moment of wij(22) in section A.5.
# This is the sample moment of h(z)^2 * t(mu(22)_i) * mu(22)_j.
# When i=j this is non-zero, and valued at the sample moment of h(z)^2 / z_j^2
# The below calculates this for all j, for each measurement and averages.
# When h(z)^2 = 0, then the value is 1, presumably because the limit of the function approaching the boundary is 1
h2onz2_mean <- function(p, n, z, h, indh){
	h4s=matrix(0,n,p)
	h4m=matrix(0,1,p)
	for (j in 1:p)
	{
		for (i in 1:n)
		{
			if (h[i] > 0){h4s[i,j]=(h[i]^2)/z[i,j]^2}
			else if (indh[i]==j){h4s[i,j]=1}
		  # for Score1 estimators the indh gives the component that is smallest (or 0 if acut constraint hit)
      # so h4s[i,j] above is set to 1 whenever h is zero for that row and j is the first component that is zero.
		}
		h4m[j]=mean(h4s[,j])
	}
	return(h4m)
}
