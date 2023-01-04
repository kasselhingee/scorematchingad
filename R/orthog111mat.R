# A helper function for creating test Astar matrices
orthogmatwith111vec <- function(){
  Q <- mclust::randomOrthogonalMatrix(3, 3)
  Rmat <- rotationmatrix(rep(1, 3), Q[, 3])
  return(Rmat %*% Q)
}
