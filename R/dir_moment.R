#' @title Method of moments estimator for Dirichlet distribution
#' @param Y Compositional data (n by p matrix)
#' @param w A vector of weights to apply to the measurements in `Y`. The length of `w` must equal the number of rows of `Y`.
#' @examples
#' Y <- MCMCpack::rdirichlet(100, c(0.1, 0.2, 1.5))
#' scorecompdir::dir_moment(Y)
#' warning("who is the author for this function?")
#' @export
dir_moment <- function(Y, w = rep(1, nrow(Y))) {
# Method of Moments estimates of the Dirichlet Distribution
temp <- dim(Y); n <- temp[1]; m <- temp[2]
# X <- cbind(X,matrix(1-apply(X,1,sum)))
cmeans <- apply(Y,2, weighted.mean, w = w)
c1sqmean <- weighted.mean(Y[,1]^2, w = w)
mom <- cmeans*(cmeans[1]-c1sqmean)/(c1sqmean - ((cmeans[1])^2))
return(matrix(mom))
}


