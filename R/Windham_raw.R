#' @title Windham transform matrix for a given parameter vector
#' @description Generates the correction matrix \eqn{\tau_c(\theta) = \tau_c \theta} for models
#' with density proportional to
#' \eqn{\exp(a(\theta) t(u))}
#' with \eqn{t(u)} a vector of sufficient statistics for a measurement \eqn{u}.
#' and \eqn{a} is *linear* function.
#' The linear assumption means that \eqn{\tau_c(\theta)} is a simple matrix operation.
#'
#' Assumes that `estimator` takes an argument paramvec if any parameter elements are fixed.
#' @param cW A vector of robustness tuning constants - the parameter vector is multiplied by these when computing the log-density of each observation for the Windham weights. For the PPI model, generate `cW` easily using [ppi_cW()] and [ppi_cW_auto()].
#' @param fpcontrol A named list of control arguments to pass to `FixedPoint::FixedPoint()` for finding the robust estimate.
#' @param ... Arguments passed to `estimator`.
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


  # extract start vector from a paramvec and paramvec_start
  if (!is.null(paramvec_start)){
    if (!is.null(extraargs$paramvec)){
      starttheta <- t_us2s(extraargs$paramvec, paramvec_start)
    } else {
      starttheta <- paramvec_start
    }
  } else { #use estimator for the start values
    starttheta <- assessment$est
  }

  # calculate isfixed
  if (!is.null(extraargs$paramvec)){isfixed = t_u2i(extraargs$paramvec)}
  else {isfixed <- rep(FALSE, length(starttheta))}

  # cW checks
  stopifnot(length(cW) == length(starttheta))
  stopifnot(is.numeric(cW))
  if (any((cW * starttheta)[isfixed] != 0)){stop("Elements of cW corresponding to fixed non-zero parameters should be zero")}
 
  # Weight correction preparation 
  originalcorrectionmethod = FALSE # use the WindhamCorrection() instead of Scealy draft method
  if (originalcorrectionmethod){
   if (length(cW) > 1){ if (var(cW[cW > 1E-10]) > (1E-10)^2){ #this check because I'm not sure what the original correction method is in the presence of a different tuning constants per value
     stop("Non-zero cW values vary, which is not supported by 'Original' Windham correction")
   }}
   inWW <- (cW > 1E-10)
   cWav <- mean(cW[cW > 1E-10]) #note that cW ~~ inWW * cWav
   thetaadjuster <- WindhamCorrection_original
  } else {
    tauc <- WindhamCorrection(cW)
    taucinv <- solve(tauc)
    cWav <- NULL  #not relevant to this correction method
    thetaadjuster <- function(newtheta, previoustheta = NULL, cW = NULL, cWav = NULL){taucinv %*% newtheta}
  }

  # functions for adding paramvec_start or otherwise to estimator arguments
  if (!assessment$paramvec_start){
    additionalargsbuilder <- function(extraargs = list(), paramvec_start = NULL){
      return(extraargs) #paramvec passed as part of extraargs
    }
  } 
  if (assessment$paramvec_start){
    additionalargsbuilder <- function(extraargs = list(), paramvec_start = NULL){
      extraargs$paramvec_start <- paramvec_start #overwrites or adds a new element to the argument list
      return(extraargs) #paramvec passed as part of extraargs
    }
  }

  # define the function that extracts the estimated parameter value
  getparamfun <- extract_paramvec_fun(assessment$estlocation)
 
  ############# 
  # build the function that takes a theta and returns a new theta, depending on assessment results
  #############
  
  fpiterator <- function(fitted){
      stopifnot(length(fitted) == sum(!isfixed))
      fulltheta <- t_sfi2u(fitted, starttheta, isfixed) #including fitted and non-fitted parameter elements
      previous <- fulltheta
      weight_vec <- Windham_weights(ldenfun = ldenfun, Y = Y,
                 theta = fulltheta, cW)

      #calculate estimate:
      args = c(list(Y = Y, w = weight_vec), additionalargsbuilder(extraargs, fulltheta)) #paramvec passed
      estobj <- do.call(estimator, args = args)
      estparamvec <- getparamfun(estobj) #extract result
      #### adjust the estimates (Step 4 in Notes5.pdf)
      estparamvec <- thetaadjuster(estparamvec, previous, cW, cWav) #for WindhamCorrections() only estparamvec is used
      fitted <- t_si2f(estparamvec, isfixed)
      return(fitted)
  }

  rlang::warn("Using the FixedPoint package - should investigate alternatives",
              .frequency = "once",
              .frequency_id = "FixedPoint_package")
  est <- fp(Function = fpiterator, Inputs = starttheta[!isfixed],
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

WindhamCorrection <- function(cW){
  weightthetamat <- diag(cW, nrow = length(cW)) #matrix that converts theta to the new theta*cW based on inclusion/exclusion  #klh: the extra argument nrow = length(cW) forces diag() to use the cW values on the diagonal, rather than treat them as the size of the matrix desired - useful when cW is legitimately length 1
  tauc <- weightthetamat + diag(1,nrow = length(cW))
  return(tauc)
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

# returns adjustment of the new estimate according to Scealy's first draft
WindhamCorrection_original <- function(newtheta, previoustheta, cW, cWav){ #cW is a vector, cWav is the average of the non-zero elements of cW
  # generate the tuning constants dbeta, dA
  dtheta <- theta * (cW - cWav) # = theta * cWav * (inWW - 1) = theta * cWav * -1 * !inWW = -cW * theta * !inWW
  return((newtheta - dtheta)/(cWav + 1))
}


