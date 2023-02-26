#' @title Build an object specifying the manifold-transform pair.
#' @description Generate an object used to specify a manifold and transformation to the manifold. This specific to each type of data.
#' @param name Name of the manifold-transform pair. See details.
#' @return An `Rcpp::XPtr` object with attribute `name`.
#' @details
#' Available pairs are
#'  + "sphere" for square-root transformation from the simplex to the positive orthant of the sphere
#'  + "simplex" for the simplex without any transformation.
#'  + "Ralr" for the additive log-ratio transformation from the simplex to Euclidean space, using the final component of vectors in the denominator of the ratio.
#'  + "Snative" for the sphere without any transformation
#' @examples
#' manifoldtransform("alr", "Euc")
#' @export
manifoldtransform <- function(start, tran, man){
  stopifnot(c(start, tran, man) %in% mantrancombos)
  out <- list(tran = new(mantranmodule$transform_ad, tran),
       man = new(mantranmodule$man_ad, man))
  return(out)
}

mantranmodule <- Rcpp::Module("manifolds", PACKAGE="scorecompdir")

mantrancombos <- matrix(c(
  "sim", "sqrt" , "sph",
  "sim", "identity", "sim",
  "sim", "alr", "Euc",
  "sim", "clr", "Hn111",
  "sph", "identity", "sph"
), byrow = TRUE, ncol = 3)
colnames(mantrancombos) <- c("start", "tran", "man")
mantrancombos <- split(mantrancombos, 1:nrow(mantrancombos))

