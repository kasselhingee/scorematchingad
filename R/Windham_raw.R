#' @title Windham transform matrix for a given parameter vector
#' @description Generates the correction matrix \eqn{\tau_c(\theta) = \tau_c \theta} for models
#' with density proportional to
#' \eqn{\exp(a(\theta) t(u))}
#' with \eqn{t(u)} a vector of sufficient statistics for a measurement \eqn{u}.
#' and \eqn{a} is *linear* function.
#' The linear assumption means that \eqn{\tau_c(\theta)} is a simple matrix operation.
#'
#' Assumes that `estimator` takes an argument paramvec if any parameter elements are fixed.
#' @export
WindhamRobust <- function(Y, estimator, ldenfun, cW, ..., fpcontrol = NULL, paramvec_start = NULL){#... earlier so that fpcontrol and paramvec_start can only be passed by being named
  extraargs <- list(...)
  ellipsis::check_dots_used()
  ### build the estimatorfun from estimator to match windham_raw() [messy but will work for testing now, then can clean later]
  # assuming estimator has arguments: Y, paramvec, w, and optionally paramvec_start.
  # and assume that the return vector can be extracted using `extract_paramvec()` and this it is the full model parameter vector, including the fixed elements.
  # also pass in ...
  ##### there are probably ways to do this where all the different types of estimatorfun are created seperately, rather than one estimator fun that does them all

  estargs <- c(list(Y = Y), extraargs)
  estargs$paramvec_start <- paramvec_start #adding this slot this way so that it is omitted if NULL
  assessment <- do.call(test_estimator2, c(list(estimator = estimator), estargs))
  estimatorfun <- function(Y, starttheta, isfixed, w){
    paramvec <- starttheta
    paramvec[!isfixed] <- NA
    args = c(list(
      Y = Y,
      w = w
      ),
      extraargs)
    if ("paramvec" %in% formalArgs(estimator)){
      args$paramvec <- paramvec
    }
    if ("paramvec_start" %in% formalArgs(estimator)){
      paramvec_start <- starttheta
      args$paramvec_start <- paramvec_start #overwrites or adds a new element to the argument list
    }
    estobj <- do.call(estimator, args = args)

    #extract estimated vector
    estparamvec <- extract_paramvec(estobj)
    return(estparamvec)
  }



  # extract start vector from a paramvec and paramvec_start
  if (!is.null(paramvec_start)){
    if (!is.null(extraargs$paramvec)){
      starttheta <- t_us2s(extraargs$paramvec, paramvec_start)
    } else {
      starttheta <- paramvec_start
    }
  } else { #use estimator for the start values
    estobj <- do.call(estimator, args = c(list(Y = Y), extraargs))
    starttheta <- extract_paramvec(estobj)
  }

  # calculate isfixed
  if (!is.null(extraargs$paramvec)){isfixed = t_u2i(extraargs$paramvec)}
  else {isfixed <- rep(FALSE, length(starttheta))}

  # check estimator fun
  chck <- estimatorfun(Y, starttheta, isfixed, w = rep(1, nrow(Y)))

  out <- windham_raw(
    prop = Y,
    cW = cW,
    ldenfun = ldenfun,
    estimatorfun = estimatorfun,
    starttheta = starttheta,
    isfixed = isfixed,
    originalcorrectionmethod = FALSE,
    fpcontrol = fpcontrol
  )
  return(out)
}

#' @param cW A vector of robustness tuning constants - the parameter vector is multiplied by these when computing the log-density of each observation for the Windham weights. For the PPI model, generate `cW` easily using [ppi_cW()] and [ppi_cW_auto()].
WindhamCorrection <- function(cW){
  weightthetamat <- diag(cW, nrow = length(cW)) #matrix that converts theta to the new theta*cW based on inclusion/exclusion  #klh: the extra argument nrow = length(cW) forces diag() to use the cW values on the diagonal, rather than treat them as the size of the matrix desired - useful when cW is legitimately length 1
  tauc <- weightthetamat + diag(1,nrow = length(cW))
  return(tauc)
}

#' @title Windham weights for a given parameter vector
#' @description Generates the weights for each measurement.
#' @param cW A robustness tuning constant. One value per element of theta
#' @param ldenfun A (possibly improper) log density function taking two arguments, `Y` and `theta`.
#' @param theta Parameters for the model
#' TRUE if the parameter is used in the Windham weights.
#' FALSE if the parameter is set to zero in the Windham weights.
#' @param sample A matrix of measurements. Each row a measurement.
#' @details
#' The elements of theta will be multiplied by cW for calculating the weights
WindhamWeights <- function(ldenfun, Y, theta, cW){
  if (is.null(ldenfun)){stop("ldenfun is NULL")}
  stopifnot(length(cW) == length(theta))
  stopifnot(is.numeric(cW))
  thetaforweights <- cW * theta #the elements of theta with FALSE inWW will be set to zero
  weights <- exp(ldenfun(Y = Y, theta = thetaforweights))
  weights=nrow(Y)*(weights/sum(weights))
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
#' @param estimatorfun A function f(Y, starttheta, isfixed, w, ...) that generates parameter estimates given
#' Y = a sample,
#' starthteta = a vector of theta values to start the iteration
#' isfixed = a vecotr of booleans. FALSE means that element of theta is estimated,
#' w = a vector of weights
#' TRUE means that element of theta is fixed at the value in `theta`.
#' The result of `estimator` must a numeric vector the same length as `starttheta`.
test_estimator <- function(estimator, Y, starttheta, isfixed, w){
  if (is.null(estimator)){stop("estimator is NULL")}
  if(!all(names(formals(estimator)[1:4]) == c("Y",  "starttheta", "isfixed", "w"))){
    stop("First four arguments of estimator must be called: Y, starttheta, isfixed, and w.")
  }

  if (is.null(w)){
    w <- runif(nrow(Y), 0, 1)
    w <- w / sum(w)
  }

  newtheta <- estimator(Y, starttheta, isfixed, w)
  if (!isTRUE(class(newtheta) == "numeric")){stop("Estimator must return a numeric value")}
  if (!isTRUE(is.vector(newtheta))){stop("Estimator must return a vector")}
  if (!isTRUE(length(newtheta) == length(starttheta))){stop("Estimator must return a vector of the same length as the input parameter vector")}
  if (any(abs(newtheta[isfixed] - starttheta[isfixed]) > sqrt(.Machine$double.eps))){
    stop("The fixed elements of theta are altered by estimator.")
  }
  invisible(NULL)
}


#' @param fpcontrol A named list of control arguments to pass to `FixedPoint::FixedPoint()` for finding the robust estimate.
#' @param ... Arguments passed to `estimator`.
#' @export
windham_raw <- function(prop, cW, ldenfun, estimatorfun, starttheta, isfixed, originalcorrectionmethod = TRUE,
   fpcontrol = NULL,
    ...)

{
  
  # cW checks
  stopifnot(length(cW) == length(starttheta))
  stopifnot(is.numeric(cW))
  if (any((cW * starttheta)[isfixed] != 0)){stop("Elements of cW corresponding to fixed non-zero parameters should be zero")}

  test_estimator(estimatorfun, prop[1:min(50, nrow(prop)), , drop = FALSE], starttheta, isfixed, w = NULL)

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
   weight_vec <- WindhamWeights(ldenfun = ldenfun, Y = prop,
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
   weight_vec <- WindhamWeights(ldenfun = ldenfun, Y = prop,
                 theta = theta, cW = cW)
   if (length(cW) > 1){ if (var(cW[cW > 1E-10]) > (1E-10)^2){ #this check because I'm not sure what the original correction method is in the presence of a different tuning constants per value
     stop("Non-zero cW values vary, which is not supported by 'Original' Windham correction")
   }}
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
