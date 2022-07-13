# A helper function for creating test Astar matrices
orthogmatwith111vec <- function(){
  Q <- mclust::randomOrthogonalMatrix(3, 3)
  Rmat <- Directional::rotation(Q[,3], rep(1, 3) / sqrt(3))
  return(Rmat %*% Q)
}
