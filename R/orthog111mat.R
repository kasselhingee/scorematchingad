#' A helper function for creating test Astar matrices
orthogmatwith111vec <- function(){
  Q <- mclust::randomOrthogonalMatrix(3, 3)
  rotmat1 <- vec2northpole(Q[,3])
  rotmat2 <- vec2northpole(rep(1, 3) / sqrt(3))
  return(rotmat2 %*% solve(rotmat1) %*% Q)
}
