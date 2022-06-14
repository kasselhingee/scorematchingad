#' @title Windham transform matrix for a given parameter vector
#' @description Generates the correction matrix \eqn{\tau_c(\theta) = \tau_c \theta} for models
#' with density proportional to
#' \eqn{\exp(a(\theta) t(u))}
#' with \eqn{t(u)} a vector of sufficient statistics for a measurement \eqn{u}.
#' and \eqn{a} is *linear* function.
#' The linear assumption means that \eqn{\tau_c(\theta)} is a simple matrix operation.
#' @param cW A robustness tuning constant.
#' @param inWW Boolean vector.
#' TRUE if the parameter is used in the Windham weights.
#' FALSE if the parameter is set to zero in the Windham weights.
WindhamCorrection <- function(cW, inWW){
  weightthetamat <- diag(inWW) #matrix that converts theta to the new theta based on inclusion/exclusion
  tauc <- cW * weightthetamat + diag(1,nrow = length(inWW))
  # tauc <- (2*cW + 1) * diag(1,nrow = length(inWW))  - cW * weightthetamat
  return(tauc)
}

#' @title Windham weights for a given parameter vector
#' @description Generates the weights for each measurement.
#' @param cW A robustness tuning constant.
#' @param ldenfun A (possibly improper) density function taking two arguments, `sample` and `theta`.
#' @param theta Parameters for the model
#' @param inWW Boolean vector.
#' TRUE if the parameter is used in the Windham weights.
#' FALSE if the parameter is set to zero in the Windham weights.
#' @param sample A matrix of measurements. Each row a measurement.
#' @details
#' The elements of theta with FALSE inWW will be set to zero for calculating the weights
WindhamWeights <- function(ldenfun, sample, theta, cW, inWW){
  thetaforweights <- inWW * theta #the elements of theta with FALSE inWW will be set to zero
  weights <- exp(ldenfun(sample, thetaforweights) * cW)
  weights=nrow(sample)*(weights/sum(weights))
  return(weights)
}

#' @title Windham weights and corrections for a given parameter vector
#' @description Generates the weights for each measurement and the correction transform \eqn{\tau_c^{-1}}
#' This function assumes that the model has density proportional to
#' \eqn{\exp(a(\theta) t(u))}
#' with \eqn{t(u)} a vector of sufficient statistics for a measurement \eqn{u}.
#' and \eqn{a} is *linear* function.
#' The linear assumption means that \eqn{\tau_c^{-1}} is a simple matrix operation.
#' @param cW A robustness tuning constant.
#' @param density A density function taking two arguments, `sample` and `theta`.
#' @param theta Parameters for the model
#' @param inWW Boolean vector.
#' TRUE if the parameter is used in the Windham weights.
#' FALSE if the parameter is set to zero in the Windham weights.
#' @param sample A matrix of measurements. Each row a measurement.



#' @title Robust score matching estimates for the generalised-gamma form of the PPI model
#' @description Uses Windham weights after alr transform to estimate parameters for the PPI model with b_L=0.
#' The final (pth) component of beta0 is not estimated.
#' @param prop compositional data (n by p matrix)
#' @param cW the robustness tuning constant c
#' @param starttheta
#' @param fixedtheta Boolean vector.
#' TRUE if parameter is fixed at the starting value.
#' FALSE if the parameter is to be estimated.
#' @param inWW Boolean vector.
#' TRUE if the parameter is used in the Windham weights.
#' FALSE if the parameter is set to zero in the Windham weights.
#' @param estimator A function that generates parameter estimates given
#' a sample,
#' a vector of weights
#' a vector of theta values (where NA are estimated, non-NA are fixed - as in buildsmotape)
#' a vector of starting values for the estimator.
#' @export
windham_raw=function(prop,cW, estimator)
{

  stop("Function empty")

}
