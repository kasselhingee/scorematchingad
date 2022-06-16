#' @title Robust score matching estimates for the generalised-gamma form of the PPI model
#' @description Uses Windham weights after alr transform to estimate parameters for the PPI model with b_L=0.
#' The final (pth) component of beta0 is not estimated.
#' @param prop compositional data (n by p matrix)
#' @param cW the robustness tuning constant c
#' @param ALs_est initial values of A_L parameter matrix
#' @param bL_est value of b_L (this should be set to zero since not estimated)
#' @param beta0_est initial values of beta (beta[p]=0 and not estimated)
#' @param ind_weightA Roughly: FALSE for the dimensions which have negative beta elements, TRUE for the dimensions with positive beta elements.
#' The pth element is not included because it is assumed that beta[p] is positive.
#' @param originalcorrectionmethod A development argument.
#' TRUE uses Janice's bias correction, FALSE uses Kassel's calculation of \eqn{tau_c} in `WindhamCorrection()`.
#' The two methods are probably equivalent, but Kassel hasn't been able to see how they are equivalent yet.
#' @export
windham_diff=function(prop,cW,ALs_est,bL_est,beta0_est, ind_weightA, originalcorrectionmethod = TRUE)
{
  stopifnot(isSymmetric(ALs_est))
  p <- ncol(prop)
  n <- nrow(prop)
  theta <- toPPIparamvec(ALs_est, bL_est, beta0_est)
	sp=p-1
	stopifnot(length(ind_weightA) == sp)
	stopifnot(all(bL_est == 0))

	ppildenfun <- function(sample, theta){
	  ppiparmats <- fromPPIparamvec(theta)
	  logden <- qldppi(sample, ppiparmats$beta, ppiparmats$ALs, ppiparmats$bL)
	  return(logden)
	}

	#indicator for each parameter of the full ppi model
	ALs_ww <- matrix(0, p-1, p-1)
	ALs_ww[!ind_weightA, !ind_weightA] <- 1
	inWW <- ppi_cppad_thetaprocessor(p, AL = ALs_ww, bL = FALSE, beta = FALSE)

	if (!originalcorrectionmethod){
  	tauc <- WindhamCorrection(cW, inWW)
  	taucinv <- solve(tauc)
	}

	stop1=0
	while(stop1==0)
	{
    # create the vector of weights
	  weight_vec <- WindhamWeights(ldenfun = ppildenfun, sample = prop,
	                 theta = theta, cW, inWW)

    if (originalcorrectionmethod){
      # generate the tuning constants dbeta, dA
      dtheta <- -cW * theta * (!inWW)
    }

    previous <- theta
    #calculate scoring estimate:
    estimator=estimatorlog_weight(prop,beta0_est[p],weight_vec)$ppi
    theta=estimator

    ### correct estimates (Step 4 in Notes5.pdf)
    if (originalcorrectionmethod){
      theta <- (theta - dtheta)/(cW+1)
    } else {
      theta <- taucinv %*% theta
    }

    #check if beta0_est has converged using the first element of the beta
    betaind <- ppithetalength(p) - p + 1
    if ( is.na(theta[betaind])) {stop("First element of beta has become NA - cannot continue")}
    if ( (abs(theta[betaind]-previous[betaind])) < 0.000001 ) {stop1=1}


	}

    estmats <- fromPPIparamvec(theta)

	return(list(ALs_est=estmats$ALs,beta0_est=estmats$beta,weight_vec=weight_vec,estimate5=theta))


}

ppiparforww <- function(beta0, ALs, bL, incinAL){
ALs_ww <- matrix(0, nrow(ALs), ncol(ALs))
ALs_ww[incinAL, incinAL] <- ALs[incinAL, incinAL]
return(list(beta0 = 0 * beta0, ALs = ALs_ww, bL = 0 * bL))
}
