#' @title Convert a PPI Parameter Vector to AL, bL and beta
#' @description
#' The parameters for the PPI model are the symmetric matrix \eqn{A_L},
#' the vector \eqn{b_L}, and the vector \eqn{\beta}.
#' However computationally, it is necessary to convert these into a single vector, which is typically performed using [`ppi_paramvec()`].
#' The function [`ppi_parammats()`] is roughly the reverse operation: 
#' any vector created by [`ppi_paramvec()`] is converted into \eqn{A_L}, \eqn{b_L} and \eqn{\beta}.
#' @param paramvec A PPI parameter vector, typically created by [`ppi_paramvec()`] or as an output of [`ppi()`].
#' @return A named list of \eqn{A_L}, \eqn{b_L}, and \eqn{\beta}.
#' @examples
#' vec <- ppi_paramvec(AL = rsymmetricmatrix(4), beta = c(-0.8, -0.7, 0))
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

