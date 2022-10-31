#' @title Generate default vector of Windham exponents for the PPI model
#' @description For the PPI model, a good default is for Windham exponents for every element of the \eqn{A_L}{A_L} matrix that corresponds to a component with concentrations near zero.
#' @param cW The value of the Windham exponent, the same for each element corresponding to the component.
#' @param ... Values of `TRUE` or `FALSE` in the same order of the components specifying if each component has a concentration near zero.
#' @return A vector of the same length as the parameter vector of the PPI model. Elements of \eqn{A_L}{A_L} will have a value of `cW` if both their row and column component has probability mass concentrated near zero. Similarly, elements of \eqn{b_L} will have a value of `cW` if their row corresponds to a component that has a concentration near zero. All other elements are zero.
#' @details
#' The final component is not included in \eqn{A_L} or \eqn{b_L} and so this function assumes that all parameters for the final component have a Windham exponent of `0`.
#' @export
ppi_cW <- function(cW, ...){
  bools <- unlist(list(...))
  bools <- as.logical(bools)
  if (any(is.na(bools))){stop("An argument after 'cW' is neither TRUE nor FALSE")}
  #if (tail(bools, 1)){stop("The final component must be set to FALSE")}

  p <- length(bools) 
  ALs_ww <- matrix(0, p-1, p-1)
  ALs_ww[bools[-p], bools[-p]] <- 1
  inWW <- ppi_paramvec(p, AL = ALs_ww, bL = bools[-p], beta = FALSE)

  cW <- cW * inWW
  return(cW)
}

#' @param Y A matrix of observations
#' @describeIn ppi_cW Automatically detects concentrations near zero by fitting a PPI distribution with \eqn{A_L=0} and \eqn{b_L=0} (i.e. a Dirichlet distribution) with the square root transformation and minima-based divergence weight function.
#' @examples
#' Y <- ppi_egmodel(100)$sample
#' ppi_cW_auto(0.01, Y)
#' ppi_cW(0.01, TRUE, TRUE, FALSE)
#' @export
ppi_cW_auto <- function(cW, Y){
  anest <- ppi(Y = Y, paramvec = ppi_paramvec(p = ncol(Y), AL=0, bL=0),
    trans = "sqrt", bdryweight = "minsq", acut = 0.01, method = "direct")
  betaest <- anest$est$beta
  cW <- do.call(ppi_cW, c(list(cW = cW), as.list(betaest < 0)))
  return(cW)
}


