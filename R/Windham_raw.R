#' @title Windham transform matrix for a given parameter vector
#' @description Generates the correction matrix \eqn{\tau_c(\theta) = \tau_c \theta} for models
#' with density proportional to
#' \eqn{\exp(a(\theta) t(u))}
#' with \eqn{t(u)} a vector of sufficient statistics for a measurement \eqn{u}.
#' and \eqn{a} is *linear* function.
#' The linear assumption means that \eqn{\tau_c(\theta)} is a simple matrix operation.
#' @param cW A robustness tuning constant.
#' TRUE if the parameter is used in the Windham weights.
#' FALSE if the parameter is set to zero in the Windham weights.
WindhamCorrection <- function(cW){
  weightthetamat <- diag(cW) #matrix that converts theta to the new theta*cW based on inclusion/exclusion
  tauc <- weightthetamat + diag(1,nrow = length(cW))
  return(tauc)
}

#' @title Windham weights for a given parameter vector
#' @description Generates the weights for each measurement.
#' @param cW A robustness tuning constant. One value per element of theta
#' @param ldenfun A (possibly improper) density function taking two arguments, `sample` and `theta`.
#' @param theta Parameters for the model
#' TRUE if the parameter is used in the Windham weights.
#' FALSE if the parameter is set to zero in the Windham weights.
#' @param sample A matrix of measurements. Each row a measurement.
#' @details
#' The elements of theta will be multiplied by cW for calculating the weights
WindhamWeights <- function(ldenfun, sample, theta, cW){
  stopifnot(length(cW) == length(theta))
  stopifnot(is.numeric(cW))
  thetaforweights <- cW * theta #the elements of theta with FALSE inWW will be set to zero
  weights <- exp(ldenfun(sample, thetaforweights))
  weights=nrow(sample)*(weights/sum(weights))
  return(weights)
}




#' @title Robust score matching estimates for the generalised-gamma form of the PPI model
#' @description Uses Windham weights after alr transform to estimate parameters for the PPI model with b_L=0.
#' The final (pth) component of beta0 is not estimated.
#' @param prop compositional data (n by p matrix)
#' @param cW the robustness tuning constant c for each element of the parameter vector
#' @param starttheta
#' @param isfixed Boolean vector.
#' TRUE if parameter is fixed at the starting value.
#' FALSE if the parameter is to be estimated.
#' TRUE if the parameter is used in the Windham weights.
#' FALSE if the parameter is set to zero in the Windham weights.
#' @param estimator A function f(Y, starttheta, isfixed, w, ...) that generates parameter estimates given
#' Y = a sample,
#' starthteta = a vector of theta values to start the iteration
#' isfixed = a vecotr of booleans. FALSE means that element of theta is estimated,
#' w = a vector of weights
#' TRUE means that element of theta is fixed at the value in `theta`.
#' @param fpcontrol A named list of control arguments to pass to `FixedPoint::FixedPoint()` for finding the robust estimate.
#' @param ... Arguments passed to `estimator`.
#' @export
windham_raw <- function(prop, cW, ldenfun, estimatorfun, starttheta, isfixed, originalcorrectionmethod = TRUE,
   fpcontrol = NULL,
    ...)

{
  if(!all(names(formals(estimatorfun)[1:4]) == c("Y",  "starttheta", "isfixed", "w"))){
    stop("First four arguments of estimatorfun must be called: Y, starttheta, isfixed, and w.")
  }

  if (!originalcorrectionmethod){
    tauc <- WindhamCorrection(cW)
    taucinv <- solve(tauc)
  }

  myfun <- function(unisfixed){
    stopifnot(length(unisfixed) == sum(!isfixed))
    previous <- unisfixed
    fulltheta <- starttheta
    fulltheta[!isfixed] <- unisfixed
    if (originalcorrectionmethod){
      theta <- Windham_raw_newtheta_original(prop = prop, cW = cW,
          ldenfun = ldenfun,
          estimatorfun = estimatorfun,
          theta = fulltheta,
          isfixed = isfixed)
    } else {
      theta <- Windham_raw_newtheta(prop = prop, cW = cW,
          ldenfun = ldenfun,
          estimatorfun = estimatorfun,
          theta = fulltheta,
          isfixed = isfixed,
          taucinv = taucinv)
    }
    unisfixed <- theta[!isfixed]
    return(unisfixed)
  }

  rlang::warn("Using the FixedPoint package - should investigate alternatives",
              .frequency = "once",
              .frequency_id = "FixedPoint_package")
  est <- fp(Function = myfun, Inputs = starttheta[!isfixed],
                    control = fpcontrol)
  nevals <- ncol(est$Inputs)
  #print(abs(est$Inputs[11, nevals] -
  #      est$Inputs[11, nevals - 1]))
  theta <- starttheta
  theta[!isfixed] <- est$FixedPoint

  return(list(theta = theta, 
           optim = list(FixedPoint = est$FixedPoint,
                        fpevals = est$fpevals,
                        Finish = est$Finish)))
}

#new theta using Kassel's correction
Windham_raw_newtheta <- function(prop, cW, ldenfun, estimatorfun, theta, isfixed, taucinv, ...){
   # create the vector of weights
   weight_vec <- WindhamWeights(ldenfun = ldenfun, sample = prop,
                 theta = theta, cW)

   #calculate estimate:
# a sample,
# a vector of weights
# a vector of theta values (where NA are estimated, non-NA are fixed - as in buildsmotape)
# a vector of starting values for the estimator.
   theta <- estimatorfun(Y = prop, starttheta = theta, isfixed = isfixed, 
                         w = weight_vec, ...)

   ### correct estimates (Step 4 in Notes5.pdf)
   theta <- taucinv %*% theta
   return(theta)
}

#new theta using original correction
Windham_raw_newtheta_original <- function(prop, cW, ldenfun, estimatorfun, theta, isfixed, ...){
   # create the vector of weights
   weight_vec <- WindhamWeights(ldenfun = ldenfun, sample = prop,
                 theta = theta, cW = cW)
   stopifnot(var(cW[cW > 1E-10]) < (1E-10)^2) #this check because I'm not sure what the original correction method is in the presence of a different tuning constants per value
   inWW <- (cW > 1E-10)
   cW <- mean(cW[cW > 1E-10])

   # generate the tuning constants dbeta, dA
   dtheta <- -cW * theta * (!inWW)
   #calculate estimate:
   theta <- estimatorfun(Y = prop, starttheta = theta, isfixed = isfixed, 
                         w = weight_vec, ...)

   ### correct estimates (Step 4 in Notes5.pdf)
   theta <- (theta - dtheta)/(cW+1)
   return(theta)
}
