#' @title Generate default vector of Windham exponents for the PPI model
#' @description For the PPI model, a good default is for Windham exponents for every element of the \eqn{A_L}{A_L} matrix that corresponds to a component with concentrations near zero.
#' @param cW The value of the Windham exponent, the same for each element corresponding to the component.
#' @param ... Values of `TRUE` or `FALSE` in the same order of the component. The final components must be included, and must be `FALSE`.
#' @return A vector of the same length as the parameter vector of the PPI model. Elements of \eqn{A_L}{A_L} will have a value of `cW` if both their row and column component has probability mass concentrated near zero. All other elements are zero.
#' @details
#' The final component is not included in \eqn{A_L} and so this function assumes that all parameters for the final component have a Windham exponent of `0`.
#' @export
ppi_cW <- function(cW, ...){
  bools <- unlist(list(...))
  bools <- as.logical(bools)
  if (any(is.na(bools))){stop("An argument after 'cW' is neither TRUE nor FALSE")}
  if (tail(bools, 1)){stop("The final component must be set to FALSE")}

  p <- length(bools) 
  ALs_ww <- matrix(0, p-1, p-1)
  ALs_ww[bools, bools] <- 1
  inWW <- ppi_paramvec(p, AL = ALs_ww, bL = FALSE, beta = FALSE)

  cW <- cW * inWW
  return(cW)
}

#' @param Y A matrix of observations
#' @describeIn ppi_cW Automatically detects concentrations near zero by fitting a Dirichlet distribution via method of moments `dir_moment`.  
#' @examples
#' Y <- MCMCpack::rdirichlet(100, c(0.1, 0.2, 1.5))
#' ppi_cW_auto(0.01, Y)
#' ppi_cW(0.01, TRUE, TRUE, FALSE)
#' @export
ppi_cW_auto <- function(cW, Y){
  betaest <- cdabyppi::dir_moment(Y)-1
  return(ppi_cW(cW, betaest < 0))
}


