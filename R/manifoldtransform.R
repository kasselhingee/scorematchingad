#' @describeIn buildsmdtape Build an object specifying the transformation from the natural domain of log-likelihood to another manifold.
#' @param start The starting manifold. Used for checking that `tran` and `man` match.
#' @param tran The name of a transformation. Available transformations are
#'  + ``srqt''
#'  + ``alr''
#'  + ``clr''
#'  + ``none'' (or ``identity'')
#' @param end The name of the manifold that `tran` maps `start` to. Available manifolds are:
#'  + ``sph'' unit sphere
#'  + ``Hn111'' hyperplane normal to the vector \eqn{1, 1, 1, 1, ...}
#'  + ``sim'' simplex
#'  + ``Euc'' Euclidean space
#' @return `manifoldtransform()` returns a named list with:
#'  + `tran` A object of type `Rcpp_transform_ad` representing the transform
#'  + `man` A object of type `Rcpp_man_ad` representing the end manifold
#' @details
#' Only some combinations of `start`, `tran` and `end` are available because `tran` must map between `start` and `end`.
#' These combinations of `start`-`tran`-`end` are currently available:
#'
#' ```{r mantrans, results = "asis", echo = FALSE}
#' cat(paste(" +", mantrancombos), sep = "\n")
#' ```
#' @examples
#' manifoldtransform("sim", "alr", "Euc")
#' @export
manifoldtransform <- function(start, tran = "identity", end = start){
  stopifnot(paste(start, tran, end, sep = "-") %in% mantrancombos)
  out <- list(tran = methods::new(mantranmodule$transform_ad, tran),
       man = methods::new(mantranmodule$man_ad, end))
  return(out)
}

mantranmodule <- Rcpp::Module("manifolds", PACKAGE="scorecompdir")

mantrancombos <- c(
  "sim-sqrt-sph",
  "sim-identity-sim",
  "sim-alr-Euc",
  "sim-clr-Hn111",
  "sph-identity-sph"
)

