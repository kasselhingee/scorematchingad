#' @title Windham Robustification of Point Estimators for Exponential Family Distributions
#' @description Performs a generalisation of Windham's robustifying method \insertCite{windham1995ro}{scorecompdir} for exponential models with natural parameters that are a linear function of the parameters for estimation.
#' Estimators must solve estimating equations of the form
#' \deqn{\sum_{i = 1}^n U(z_i; \theta) = 0.}
#' The estimate is found iteratively through a fixed point method as suggested by \insertCite{windham1995ro;textual}{scorecompdir}.
#'

#' @param Y A matrix of measurements. Each row is a measurement, each component is a dimension of the measurement.
#' @param estimator A function that estimates parameters from weighted observations.
#' It must have arguments `Y` that is a matrix of measurements and `w` that are weights associated with each row of `Y`. If it accepts arguments `paramvec` or `paramvec_start` then these will be used to specify fixed elements of the parameter vector and the starting guess of the parameter vector, respectively. The estimated parameter vector, including any fixed elements, must be the returned object, or the first element of a returned list, or as the `paramvec` slot within the `est` slot of the returned object.
#' @param ldenfun A function that returns a vector of values propotional to the log-density for a matrix of observations `Y` and parameter vector `theta`.
#' @param ... Arguments passed to `estimator`.
#' @param paramvec_start
#' Initially used to check the function `estimator`. If `estimator` accepts a `paramvec_start`, then the current estimate of the parameter vector is passed as `paramvec_start` to `estimator` in each iteration.
#' @param cW A vector of robustness tuning constants. When computing the weight for an observation the parameter vector is multiplied element-wise with `cW`. For the PPI model, generate `cW` easily using [ppi_cW()] and [ppi_cW_auto()].
#' @param fpcontrol A named list of control arguments to pass to [FixedPoint::FixedPoint()] for the fixed point iteration.


#' @details
#' For any family of models with density \eqn{f(z; \theta)}, Windham's method finds the parameter set \eqn{\hat\theta} such that the estimator applied to observations weighted by \eqn{f(z; \hat\theta)^c} returns an estimate that matches the theoretical effect of weighting the full population of the model.
#' When \eqn{f} is proportional to \eqn{\exp(\eta(\theta) \cdot T(z))} and \eqn{\eta(\theta)} is linear, then these weights are equivalent to \eqn{f(z; c\hat\theta)} and the theoretical effect of the weighting on the full population is to scale the parameter vector \eqn{\theta} by \eqn{1+c}.
#' The function `Windham()` assumes that this is the case, and allows a generalisation where \eqn{c} is a vector so the weight for an observation \eqn{z} is \deqn{f(z; c \circ \theta),} where \eqn{f} is the proposed model density function,  \eqn{\theta} is the parameter vector, \eqn{c} is a vector of tuning constants, and \eqn{\circ} is the element-wise product (Hadamard product).
#'
#' The solution is found iteratively \insertCite{windham1995ro}{scorecompdir}. 
#' Given a parameter set \eqn{\theta_n}, `Windham()` first computes weights \eqn{f(z; c \circ \theta_n)} for each observation \eqn{z}.
#' Then, a new parameter set \eqn{\tilde{\theta}_{n+1}} is estimated by `estimator` with the computed weights.
#' This new parameter set is element-wise-multiplied by the (element-wise) reciprical of \eqn{1+c} to obtain an adjusted parameter set \eqn{\theta_{n+1}}.
#' The estimate returned by `Windham()` is the parameter set \eqn{\hat{\theta}} such that \eqn{\theta_n \approx \theta_{n+1}}.
#'
#' An exponential model with a base rate may be used with `Windham()` so long as the base rate is omitted from `ldenfun` (i.e. not used for weighting).
#' @family Windham functions
#' @return
#' A list:
#' * `theta` the estimated parameter vector
#' * `optim` information about the fixed point iterations and opimisation process. Including a slot `finalweights` for the weights in the final iteration.
#' @export
Windham <- function(Y, estimator, ldenfun, cW, ..., fpcontrol = NULL, paramvec_start = NULL){#... earlier so that fpcontrol and paramvec_start can only be passed by being named
  out <- Windham_raw(Y, estimator, ldenfun, cW, ..., fpcontrol = fpcontrol, paramvec_start = paramvec_start, multiplicativecorrection = TRUE)
  return(out)
}

Windham_raw <- function(Y, estimator, ldenfun, cW, ..., fpcontrol = NULL, paramvec_start = NULL, multiplicativecorrection = TRUE){#... earlier so that fpcontrol and paramvec_start can only be passed by being named
  extraargs <- list(...)
  ellipsis::check_dots_used()
  # assuming estimator has arguments: Y, paramvec, w, and optionally paramvec_start.
  # and assume that the return vector can be extracted using `extract_paramvec()` and similar and that this is the full model parameter vector, including the fixed elements (this is important for computing density).
  estargs <- c(list(Y = Y), extraargs)
  estargs$paramvec_start <- paramvec_start #adding this slot this way so that it is omitted if NULL

  #assess the passes estimator
  assessment <- do.call(Windham_assess_estimator, c(list(estimator = estimator), estargs))


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

  # Correction of parameter preparation
  # use the WindhamCorrection(), the alternative is Scealy's original additive method in the draft paper
  if (!multiplicativecorrection){
   message("Using the bias correction of Scealy et al 2023")
   if (length(cW) > 1){ if (var(cW[cW > 1E-10]) > (1E-10)^2){ #require constant cW (or zero) because I'm not sure what Scealy's correction method should be in the presence of a different tuning constants per value
     stop("Non-zero cW values vary, which is not supported by 'additive' correction of the parameter estimate.")
   }}
   inWW <- (cW > 1E-10)
   cWav <- mean(cW[cW > 1E-10]) #note that cW ~~ inWW * cWav
   thetaadjuster <- WindhamCorrection_additive
  } else {
    tauinv <- WindhamCorrection_multiplicative_tauinv(cW)
    cWav <- NULL  #not relevant to this correction method
    thetaadjuster <- function(newtheta, previoustheta = NULL, cW = NULL, cWav = NULL){tauinv %*% newtheta}
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
  theta <- starttheta
  theta[!isfixed] <- est$FixedPoint

  # get weights corresponding the final iteration
  thetaprevious <- t_sfi2u(est$Inputs[,ncol(est$Inputs)], starttheta, isfixed) #including fitted and non-fitted parameter elements
  weight_vec <- Windham_weights(ldenfun = ldenfun, Y = Y,
                                theta = thetaprevious, cW)

  return(list(theta = theta,
           optim = list(FixedPoint = est$FixedPoint,
                        fpevals = est$fpevals,
                        Finish = est$Finish,
                        finalweights = weight_vec)))
}

# @title Windham transform matrix for a given parameter vector
# @description Generates the correction matrix $\eqn{\tau_c(\theta) = \tau_c \theta}$ for models
# with density proportional to
# $\eqn{\exp(\eta(\theta) t(u))}$
# with \eqn{t(u)} a vector of sufficient statistics for a measurement \eqn{u}.
# and \eqn{\eta} is *linear* function.
# The linear assumption means that \eqn{\tau_c(\theta)} is a multiplication by the matrix diag(1 + cW).
WindhamCorrection_multiplicative_tauinv <- function(cW){
  tauinv <- diag(1/(1 + cW), nrow = length(cW)) #matrix that converts theta to the new theta*cW based on inclusion/exclusion  #klh: the extra argument nrow = length(cW) forces diag() to use the cW values on the diagonal, rather than treat them as the size of the matrix desired - useful when cW is legitimately length 1
  return(tauinv)
}

# returns adjustment of the new estimate according to Scealy's first draft
WindhamCorrection_additive <- function(newtheta, previoustheta, cW, cWav){ #cW is a vector, cWav is the average of the non-zero elements of cW
  # generate the tuning constants dbeta, dA
  dtheta <- previoustheta * (cW - cWav) # = theta * cWav * (inWW - 1) = theta * cWav * -1 * !inWW = -cW * theta * !inWW
  return((newtheta - dtheta)/(cWav + 1))
}


