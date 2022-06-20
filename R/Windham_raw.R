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
#' @param estimator A function f(Y, weights, theta, fixedtheta) that generates parameter estimates given
#' a sample Y,
#' a vector of weights
#' a vector of theta values
#' a boolean vector. FALSE means that element of theta is estimated.
#' TRUE means that element of theta is fixed at the value in `theta`.
#' @param ... Arguments passed to `estimator`.
#' @export
windham_raw <- function(prop, cW, ldenfun, estimatorfun, starttheta, fixedtheta, inWW, originalcorrectionmethod = TRUE, ...)
{
  if(!all(names(formals(estimatorfun)[1:4]) == c("Y", "weights", "starttheta", "fixedtheta"))){
    stop("First four arguments of estimatorfun must be called: Y, weights, starttheta and fixedtheta.")
  }
  #all(inWW * fixedtheta == 0) #every fixed (not fitted) parameter should not be in the Windham Weight

  if (!originalcorrectionmethod){
    tauc <- WindhamCorrection(cW, inWW)
    taucinv <- solve(tauc)
  }

  myfun <- function(unfixedtheta){
    stopifnot(length(unfixedtheta) == sum(!fixedtheta))
    previous <- unfixedtheta
    fulltheta <- starttheta
    fulltheta[!fixedtheta] <- unfixedtheta
    if (originalcorrectionmethod){
      theta <- Windham_raw_newtheta_original(prop = prop, cW = cW,
          ldenfun = ldenfun,
          estimatorfun = estimatorfun,
          theta = fulltheta,
          fixedtheta = fixedtheta,
          inWW = inWW)
    } else {
      theta <- Windham_raw_newtheta(prop = prop, cW = cW,
          ldenfun = ldenfun,
          estimatorfun = estimatorfun,
          theta = fulltheta,
          fixedtheta = fixedtheta,
          inWW = inWW,
          taucinv = taucinv)
    }
    unfixedtheta <- theta[!fixedtheta]
    return(unfixedtheta)
  }

  rlang::warn("Using the FixedPoint package - should investigate alternatives",
              .frequency = "once",
              .frequency_id = "FixedPoint_package")
  est <- FixedPoint::FixedPoint(myfun, starttheta[!fixedtheta],
                    Method = "VEA",
                    ConvergenceMetricThreshold = 1E-10)
  nevals <- ncol(est$Inputs)
  print(abs(est$Inputs[11, nevals] -
        est$Inputs[11, nevals - 1]))
  theta <- starttheta
  theta[!fixedtheta] <- est$FixedPoint

  return(list(theta = theta, 
           optim = list(FixedPoint = est$FixedPoint,
                        fpevals = est$fpevals,
                        Finish = est$Finish)))
}


#new theta using Kassel's correction
Windham_raw_newtheta <- function(prop, cW, ldenfun, estimatorfun, theta, fixedtheta, inWW, taucinv, ...){
   # create the vector of weights
   weight_vec <- WindhamWeights(ldenfun = ldenfun, sample = prop,
                 theta = theta, cW, inWW)

   #calculate estimate:
# a sample,
# a vector of weights
# a vector of theta values (where NA are estimated, non-NA are fixed - as in buildsmotape)
# a vector of starting values for the estimator.
   theta <- estimatorfun(Y = prop, weights = weight_vec,
                starttheta = theta, fixedtheta = fixedtheta, ...)

   ### correct estimates (Step 4 in Notes5.pdf)
   theta <- taucinv %*% theta
   return(theta)
}

#new theta using original correction
Windham_raw_newtheta_original <- function(prop, cW, ldenfun, estimatorfun, theta, fixedtheta, inWW, ...){
   # create the vector of weights
   weight_vec <- WindhamWeights(ldenfun = ldenfun, sample = prop,
                 theta = theta, cW, inWW)

   # generate the tuning constants dbeta, dA
   dtheta <- -cW * theta * (!inWW)
   #calculate estimate:
   theta <- estimatorfun(Y = prop, weights = weight_vec,
                starttheta = theta, fixedtheta = fixedtheta, ...)

   ### correct estimates (Step 4 in Notes5.pdf)
   theta <- (theta - dtheta)/(cW+1)
   return(theta)
}
