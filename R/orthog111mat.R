


orthogmatwith111vec <- function(){
  Q <- mclust::randomOrthogonalMatrix(3, 3)
  rotmat <- to111mat(Q[,3])
  return(rotmat %*% Q)
}

rot2d <- function(theta){
  out <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), ncol = 2, nrow = 2, byrow = TRUE)
  return(out)
}
rot3dx <- function(theta){
  out <- matrix(0, ncol = 3, nrow = 3)
  out[2:3, 2:3] <- rot2d(theta)
  out[1,1] <- 1
  return(out)
}
rot3dy <- function(theta){
  out <- matrix(0, ncol = 3, nrow = 3)
  out[c(1, 3), c(1, 3)] <- t(rot2d(theta))
  out[2,2] <- 1
  return(out)
}
rot3dz <- function(theta){
  out <- matrix(0, ncol = 3, nrow = 3)
  out[1:2, 1:2] <- rot2d(theta)
  out[3,3] <- 1
  return(out)
}

# determine angles
anglestoxaxis <- function(vec){
  ## first project on yzplane
  vecprojzyplane <- vec - sum(vec * c(1, 0, 0)) * c(1, 0, 0)
  theta1 <- atan(vecprojzyplane[3]/ vecprojzyplane[2])
  vecrotxyplane <- rot3dx(-theta1) %*% vec #use the angle to get rotate the vector onto the xyplane
  rotmat <- rot3dx(-theta1)
  stopifnot(abs(sum(vecrotxyplane^2) - 1) < 1E-10)

  #now rotate the vector onto the x-axis
  theta2 <- atan(vecrotxyplane[2]/ vecrotxyplane[1])
  xonly <- rot3dz(-theta2) %*% vecrotxyplane
  rotmat <- rot3dz(-theta2) %*% rotmat
  if (xonly[[1]] < 0){#rotate pi
    xonly <- rot3dz(pi) %*% xonly
    rotmat <- rot3dz(pi) %*% rotmat
  }
  stopifnot(abs(sum(xonly^2) - 1) < 1E-10)
  stopifnot(abs(xonly[[1]] - 1) < 1E-10)
  return(rotmat)#about the z axis
}
to111mat <- function(vec){
  rotmat <- anglestoxaxis(vec)
  rotmatOnes <- anglestoxaxis(c(1,1,1) / sqrt(3))
  out <- solve(rotmatOnes) %*% rotmat
  #check
  dir1 <- out %*% vec
  stopifnot(isTRUE(all.equal(dir1[1:2] - dir1[3], rep(0, 2))))
  return(out)
}
