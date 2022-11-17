#' @title Evaluate the Hessian and gradient offset of a Quadratic Score Matching Objective Function
#' @description
#' When the score matching objective function is quadratic then score matching objective for a single measurement can be written
#' \deqn{\tilde\psi_{f_\theta, z} = \frac{1}{2} \theta^T W(z) \theta + d(z)^T\theta.}
#' The function `quadtraticsmo_parts_approx()` evaluates \eqn{W(z)} and \eqn{d(z)} using Taylor approximation.
#' @details
#' Tapes of the Hessian and gradient offset are created using [`pTapeHessian()`] and [`pTapeGradOffset()`] respectively.
#' These tapes are then evaluated for every row of `Y`.
#' 
#' @param smotape A tape of a score matching objective function.
#' @param Y A matrix of (multivariate) observations. Each row is an observation.
quadraticsmo_parts <- function(smotape, Y){

stopifnot(testquadratictape(smotape))


}
