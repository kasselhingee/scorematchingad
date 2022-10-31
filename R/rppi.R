#' @title Generate random observations from the PPI model
#' @description Given parameters of the PPI model, generates independent samples.
#' @param n Sample size
#' @param p Dimension (number of components)
#' @param beta The \eqn{\beta}{beta} shape parameter vector
#' @param AL The \eqn{A_L} parameter matrix
#' @param bL The \eqn{b_L} parameter vector
#' @param paramvec The PPI parameter vector, created easily using [ppi_paramvec()] and also returned by [ppi()].
#' @param maxden This is the constant \eqn{log(C)} in (Scealy and Wood, 2021; Appendix A.1.1)
#' @return A matrix with `n` rows. Each row is a independent draw from the specified PPI distribution.
#' @inherit ppi sections
#' @details
#' We recommend running `rppi()` a number of times to ensure the choice of `maxden` is good. `rppi()` will error when `maxden` is too low.
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
#' AL=-0.5*solve(SigA)
#' bL=solve(SigA)%*%muL
#'
#' samp <- rppi(n,p,beta,AL,bL,4)
#' plot(ks::kde(samp[,-p]),
#'  xlim = c(0, 1), ylim = c(0, 1))
#' segments(0, 0, 0, 1)
#' segments(0, 1, 1, 0)
#' segments(1, 0, 0, 0)
#'
#' dppi(samp, AL=AL, bL=bL, beta=beta)
#' @export
rppi <- function(n, beta = NULL, AL = NULL, bL = NULL, paramvec = NULL, maxden = 4){
  # a warning if maxden is high
  if (maxden > 10){
    rlang::warn(message = paste(sprintf("'maxden' of %0.2f is higher than 10.", maxden),
                                "When rppi() requires a high 'maxden' it could mean that",
                                "PPI density is hugely different from the Dirichlet component of the density.",
                                "This could mean that the concentrations on the boundary from the Dirichlet component",
                                "will be too narrow to be represented in simulatad samples."),
                .frequency = "once",
                .frequency_id = "highmaxden")
  }

  #process inputs
  if (is.null(paramvec)){if (any(is.null(beta), is.null(AL), is.null(bL))){stop("If paramvec isn't supplied then beta, AL, and bL must be supplied")}}
  else {
    if (!all(is.null(beta), is.null(AL), is.null(bL))){stop("Providing a paramvec is incompatible with providing a beta, AL or bL.")}
    if (any(is.na(paramvec))){stop("All elements of paramvec must be non-NA")}
    parammats <- fromPPIparamvec(paramvec)
    beta <- parammats$beta
    AL <- parammats$ALs
    bL <- parammats$bL
  }
  p <- length(beta)

  maxdenin <- maxden
  # first simulate starting with a block of Dirichlet samples of size n.
  firstaccepted <- rppi_block(n,p,beta = beta,AL = AL,bL,maxden)
  maxden <- firstaccepted$maxden
  samples <- firstaccepted$accepted
  propaccepted <- max(nrow(samples) / n, 1E-3)
  # based on the number of samples that came out, simulate the remaining
  while (nrow(samples) < n){
    newsamples <- rppi_block(ceiling((n - nrow(samples)) * 1/propaccepted),p,beta,AL,bL,maxden)
    maxden <- newsamples$maxden
    samples <- rbind(samples, newsamples$accepted)
    # continue until n or more samples accepted
  }

  #maxden is the constant log(C) in Appendix A.1.3. Need to run the sampler
  #a few times to check that it is an appropriate upper bound.
  if (maxden > maxdenin){stop(sprintf("Input maxden (i.e. the log(C) maximum) was %0.2f, but sampler suggests higher than %0.2f is required.", maxdenin, maxden))}


  samples <- samples[1:n, ] # remove extra samples

  return(samples)
}


# simulates samples one at a time
rppi_singly <- function(n,p,beta,AL,bL,maxden)
{

	alpha=beta+1
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
    num <- ppi_uAstaru(matrix(Uni_nop, nrow = 1), AL, bL) - maxden #uT * AL * u + t(bL) * u - maxden
		if (num > 0){maxden=num+maxden}
		#print(maxden)
		if (u < exp(num)){samp1=rbind(samp1,Uni);coun=coun+1}
		#print(coun)
	}
	samp3=samp1[2:length(samp1[,1]),]

	return(list(samp3=samp3,maxden=maxden))
}



# try generating samples without a 'while' loop
rppi_block <- function(n,p,beta,AL,bL,maxden){
  Uni <- MCMCpack::rdirichlet(n, beta+1)
  Uni_nop <- Uni[, -p]
  nums <- ppi_uAstaru(Uni_nop, AL, bL) - maxden #uT * AL * u + t(bL) * u - maxden

  # update maxdens
  max_num <- max(nums)
  if (max_num > 0) {maxden=max_num+maxden}

  # accept some of points
  unif <- stats::runif(n,0,1)
  accepted <- Uni[unif < exp(nums), , drop = FALSE]

  return(list(accepted = accepted, maxden = maxden))
}

#' @title Improper Log-Density of the PPI Model
#' @description Compute the natural logarithm of the improper density for the PPI model for the given matrix of measurements `prop`. Rows with negative values or with a sum that is more than `1E-15` from `1` are assigned a value of `-Inf`.
#' @param Y A matrix of measurements.
#' @inheritParams rppi
#' @details The value calculated by `dppi` is
#' \deqn{z^TA_Lz + b_L^Tz + \beta^T \log(z).}
#' @export
dppi <- function(Y, AL = NULL,bL = NULL, beta = NULL, paramvec = NULL){
  #process inputs
  if (is.null(paramvec)){if (any(is.null(beta), is.null(AL), is.null(bL))){stop("If paramvec isn't supplied then beta, AL, and bL must be supplied")}}
  else {
    if (!all(is.null(beta), is.null(AL), is.null(bL))){stop("Providing a paramvec is incompatible with providing a beta, AL or bL.")}
    if (any(is.na(paramvec))){stop("All elements of paramvec must be non-NA")}
    parammats <- fromPPIparamvec(paramvec)
    beta <- parammats$beta
    AL <- parammats$ALs
    bL <- parammats$bL
  }

  p <- ncol(Y)
  if (!("matrix" %in% class(bL))){bL <- as.matrix(bL, ncol = 1)}
  stopifnot(isTRUE(ncol(bL) == 1))

  uAstaru <- ppi_uAstaru(Y[,-p, drop = FALSE], AL, bL) #result is a vector
  if (all(beta == 0)){return(as.vector(uAstaru))} #skip the computation below when beta0 is zero

  if (!("matrix" %in% class(beta))){beta <- as.matrix(beta, ncol = 1)}
  logprop <- log(Y)
  # define u^0 as 1 when u goes to -Inf
  logprop[, beta == 0] <- 0
  logdirichlet <- logprop %*% beta
  logdensity <- as.vector(uAstaru + logdirichlet)

  # set points outside the simplex to 0
  negatives <- rowSums(Y < 0) > 0
  sumisnt1 <- abs(rowSums(Y) -1) > 1E-15
  logdensity[negatives|sumisnt1] <- -Inf
  return(logdensity)
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
