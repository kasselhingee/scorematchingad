#' @title Method of moments estimator for Dirichlet distribution
#' @param X Compositional data (n by p matrix)
#' @export
dirichmom <- function(X) {
# Method of Moments estimates of the Dirichlet Distribution
temp <- dim(X); n <- temp[1]; m <- temp[2]
# X <- cbind(X,matrix(1-apply(X,1,sum)))
mom <- apply(X,2,mean)*(mean(X[,1])-mean(X[,1]^2))/(mean(X[,1]^2) - ((mean(X[,1]))^2))
return(matrix(mom))
}


