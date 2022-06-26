#' @title Generate random observations from the PPI model
#' @description Given parameters of the PPI model, generates independent samples.
#' @param n Sample size
#' @param p Dimension (number of components)
#' @param beta0 The \eqn{\beta_0}{beta0} shape parameter vector
#' @param ALs The \eqn{A_L} parameter matrix
#' @param bL The \eqn{b_L} parameter vector
#' @param maxden This is the constant \eqn{log(C)} in (Scealy and Wood, 2021; Appendix A.1.1)
#' @return A list. The first element is the sample in the form of a matrix with `n` rows and `p` columns.
#' The second element is the maxden updated based on whether the sample exceeds the input maxden.
#' @details
#' \eqn{A_L} controls the covariance between components.
#' \eqn{b_L} controls the location of the distribution within the simplex
#' \eqn{\beta_0[i]}{beta0[i]} controls the shap of the density when the ith component is close to zero.
#' If \eqn{b_L} is such that the ith component is typically far from zero, then \eqn{\beta_0[i]}{beta0[i]} will have negligible effect.
#' @examples
#' n=1000
#' p=3
#' beta0=c(-0.8, -0.8, -0.5)
#'
#' muL=matrix(0,p-1,1)
#' muL[1:sum(p,-1)]=0.12
#' aa=matrix(1/500,p-1,1)
#' D=diag(as.vector(aa))
#' SigA=D
#' SigA[1,1]=SigA[1,1]*2
#' cor=0.5
#' SigA[1,2]=cor*sqrt(SigA[1,1]*SigA[2,2])
#' SigA[2,1]=SigA[1,2]
#' ALs=-0.5*solve(SigA)
#' bL=solve(SigA)%*%muL
#'
#' samp <- rhybrid(n,p,beta0,ALs,bL,4)
#' plot(ks::kde(samp$samp3[,-p]),
#'  xlim = c(0, 1), ylim = c(0, 1))
#' segments(0, 0, 0, 1)
#' segments(0, 1, 1, 0)
#' segments(1, 0, 0, 0)
#'
#' qldppi(samp$samp3, beta0, ALs, bL)
#' @export
rhybrid <- function(n,p,beta0,ALs,bL,maxden){
  # a warning if maxden is high
  if (maxden > 10){
    rlang::warn(message = paste(sprintf("'maxden' of %0.2f is higher than 10.", maxden),
                                "When rhybrid() requires a high 'maxden' it could mean that",
                                "PPI density is hugely different from the Dirichlet component of the density.",
                                "This could mean that the concentrations on the boundary from the Dirichlet component",
                                "will be too narrow to be represented in simulatad samples."),
                .frequency = "once",
                .frequency_id = "highmaxden")
  }

  maxdenin <- maxden
  # first simulate starting with a block of Dirichlet samples of size n.
  firstaccepted <- rhybrid_block(n,p,beta0,ALs,bL,maxden)
  maxden <- firstaccepted$maxden
  samples <- firstaccepted$accepted
  propaccepted <- max(nrow(samples) / n, 1E-3)
  # based on the number of samples that came out, simulate the remaining
  while (nrow(samples) < n){
    newsamples <- rhybrid_block(ceiling((n - nrow(samples)) * 1/propaccepted),p,beta0,ALs,bL,maxden)
    maxden <- newsamples$maxden
    samples <- rbind(samples, newsamples$accepted)
    # continue until n or more samples accepted
  }

  #maxden is the constant log(C) in Appendix A.1.3. Need to run the sampler
  #a few times to check that it is an appropriate upper bound.
  if (maxden > maxdenin){stop(sprintf("Input maxden (i.e. the log(C) maximum) was %0.2f, but sampler suggests higher than %0.2f is required.", maxdenin, maxden))}


  samples <- samples[1:n, ] # remove extra samples

  return(list(samp3=samples,maxden=maxden))
}


# simulates samples one at a time
rhybrid_singly <- function(n,p,beta0,ALs,bL,maxden)
{

	alpha=beta0+1
	coun=0
	samp1=matrix(0,1,p)
	count2=0
	while (coun < n)
	{

		count2=count2+1
		Uni=MCMCpack::rdirichlet(1, alpha)
		u=stats::runif(1,0,1)
		Uni_nop <- Uni[1:(p-1)]
		tUni_nop <- t(Uni[1:(p-1)])
    num <- ppi_uAstaru(matrix(Uni_nop, nrow = 1), ALs, bL) - maxden #uT * ALs * u + t(bL) * u - maxden
		if (num > 0){maxden=num+maxden}
		#print(maxden)
		if (u < exp(num)){samp1=rbind(samp1,Uni);coun=coun+1}
		#print(coun)
	}
	samp3=samp1[2:length(samp1[,1]),]

	return(list(samp3=samp3,maxden=maxden))
}



# try generating samples without a 'while' loop
rhybrid_block <- function(n,p,beta0,ALs,bL,maxden){
  Uni <- MCMCpack::rdirichlet(n, beta0+1)
  Uni_nop <- Uni[, -p]
  nums <- ppi_uAstaru(Uni_nop, ALs, bL) - maxden #uT * ALs * u + t(bL) * u - maxden

  # update maxdens
  max_num <- max(nums)
  if (max_num > 0) {maxden=max_num+maxden}

  # accept some of points
  unif <- stats::runif(n,0,1)
  accepted <- Uni[unif < exp(nums), , drop = FALSE]

  return(list(accepted = accepted, maxden = maxden))
}

#' @describeIn rhybrid Compute the logarithm of the improper density for the PPI model for the given matrix of measurements `prop`.
#' @param `prop` A matrix of measurements.
#' @export
qldppi <- function(prop,beta0,ALs,bL){
  p <- ncol(prop)
  sp <- p - 1
  if (!("matrix" %in% class(bL))){bL <- as.matrix(bL, ncol = 1)}
  stopifnot(isTRUE(ncol(bL) == 1))
  uAstaru <- ppi_uAstaru(prop[,-p], ALs, bL) #result is a vector
  if (all(beta0 == 0)){return(as.vector(uAstaru))} #skip the computation below when beta0 is zero
  if (!("matrix" %in% class(beta0))){beta0 <- as.matrix(beta0, ncol = 1)}
  logprop <- log(prop)
  # define u^0 as 1 when u goes to -Inf
  logprop[, beta0 == 0] <- 0
  logdirichlet <- logprop %*% beta0
  return(as.vector(uAstaru + logdirichlet))
}

# below is function for uT * ALs * u + t(bL) * u
# prop_nop is the measurements WITHOUT the last column
ppi_uAstaru <- function(prop_nop, ALs, bL){
  stopifnot(ncol(prop_nop) == ncol(ALs))
  stopifnot(ncol(ALs) == nrow(bL))
  stopifnot(ncol(ALs) == nrow(ALs))
  nums <- rowSums((prop_nop %*% ALs) * prop_nop) + #uT * ALs * u
    as.vector(prop_nop %*% bL) #+ u * bL
  return(nums)
}
