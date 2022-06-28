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
  # preperation - tranforming the inputs
  p <- ncol(prop)
  theta <- toPPIparamvec(ALs_est, bL_est, beta0_est)
	sp=p-1
	stopifnot(length(ind_weightA) == sp)
	stopifnot(all(bL_est == 0))

  #indicator for each parameter of the full ppi model
  ALs_ww <- matrix(0, p-1, p-1)
  ALs_ww[!ind_weightA, !ind_weightA] <- 1
  inWW <- ppi_cppad_thetaprocessor(p, AL = ALs_ww, bL = FALSE, beta = FALSE)

  #preparing ppi specific info
  ppildenfun <- function(sample, theta){
    ppiparmats <- fromPPIparamvec(theta)
    logden <- qldppi(sample, ppiparmats$beta, ppiparmats$ALs, ppiparmats$bL)
    return(logden)
  }

  ppiestimator <- function(Y, starttheta, isfixed, w){
          estimatorlog_weight(prop = prop, betap = starttheta[length(starttheta)], weightW = w)$ppi}

  isfixed <- ppi_cppad_thetaprocessor(p, AL=FALSE, bL = TRUE, betaL = FALSE, betap = TRUE)
  est <- windham_raw(prop = prop,
               cW = cW,
               ldenfun = ppildenfun,
               estimatorfun = ppiestimator,
               starttheta = theta,
               isfixed = isfixed,
               inWW = inWW,
               originalcorrectionmethod = originalcorrectionmethod)

  estmats <- fromPPIparamvec(est$theta)

  out <- list(est = c(list(theta = est$theta), estmats),
       optim = est$optim)
                  
  return(out)
}

