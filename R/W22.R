
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
# When i not equal j then this is zero
h2onz2_mean <- function(p, n, z, h, indh){
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
	return(h4m)
}
