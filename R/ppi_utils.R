#' @noRd
#' @title To or From vector form of parameters for PPI
#' purely for testing hardcoded ppi estimators - returns canonical parameters
toPPIcannparam <- function(ALs, bL, beta, manifold = "sphere"){
  if (manifold == "sphere"){out <- ppi_paramvec(AL = ALs, bL = bL, beta = 1 + 2 * beta)}
  else if (manifold == "simplex"){out <- ppi_paramvec(AL = ALs, bL = bL, beta = beta)}
  else {stop("manifold not supported")}
  return(out)
}

#' @noRd
#' @title returns the length of the parameter vector for given dimension p
ppithetalength <- function(p){
  p + #the beta
  (p-1) + #the diagonal of AL
  (p-2) * (p-1)/2 + #the upper triangle of AL
  (p-1) #the bL
}

#' @noRd
#' @title From the length of a ppi parameter vector, get the number of components
ppiltheta2p <- function(ltheta){#ltheta is length of theta
  #ltheta = p + (p-1) + (p-2) * (p-1)/2 + (p-1)
  # = 3p - 2 + (p^2 - 3p + 2)/2
  # 0 = p^2/2 + 1.5p -1-ltheta
  # 0 = p^2 + 3p - (2+2ltheta)
  #p = (-3 pm sqrt(9 + 4(2+2ltheta)) / 2
  #p = (-3 + sqrt(9 + 8 + 8ltheta)) /2
  p <- (-3 + sqrt(9 + 8 + 8 * ltheta))/2
  return(p)
}

#' @title Convert a PPI Parameter Vector to AL, bL and beta
#' @family PPI model tools
#' @description
#' The parameters for the PPI model are the symmetric matrix \eqn{A_L},
#' the vector \eqn{b_L}, and the vector \eqn{\beta}.
#' However computationally, it is necessary to convert these into a single vector, which is typically performed using [`ppi_paramvec()`].
#' The function [`ppi_parammats()`] is roughly the reverse operation: 
#' any vector created by [`ppi_paramvec()`] is converted into \eqn{A_L}, \eqn{b_L} and \eqn{\beta}.
#' @param paramvec A PPI parameter vector, typically created by [`ppi_paramvec()`] or as an output of [`ppi()`].
#' @return A named list of \eqn{A_L}, \eqn{b_L}, and \eqn{\beta}.
#' @examples
#' vec <- ppi_paramvec(AL = rsymmetricmatrix(2), beta = c(-0.8, -0.7, 0))
#' ppi_parammats(vec)
#' @export
ppi_parammats <- function(paramvec){
  calcp <- ppiltheta2p(length(paramvec))
  p <- calcp
  AL <- tosmatrix(paramvec[1:((p-1) + (p-1)*(p-2)/2)])
  bL <- paramvec[p - 1 + ((p-2) * (p-1)/2) + 1:(p-1)]
  beta <- paramvec[(p - 1 + ((p-2) * (p-1)/2) + (p-1) + 1):length(paramvec)]
  return(list(
    AL = AL,
    bL = bL,
    beta = beta
  ))
}

