# tools for converting symmetric matrices to vectors and back

tosmatrix <- function(theta){
  # l = d + d * (d-1)/2
  # 2l = 2d + d^2 - d
  # 0 = d^2 + d - 2l
  # d = (-1 pm sqrt(1-4*(-2l)))/2
  #   = (-1 + sqrt(1+8l))/2
  d <- (-1+sqrt(8 * length(theta) + 1)) / 2
  mat <- matrix(NA, d, d)
  diag(mat) <- theta[1:d]
  mat[upper.tri(mat)] <- theta[(d + 1):length(theta)]
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  return(mat)
}

fromsmatrix <- function(mat){
  stopifnot(isSymmetric.matrix(mat))
  return(c(diag(mat), mat[upper.tri(mat)]))
}

# return the array of row/column indices, and the length
# for the upper triangle
indexcombinations <- function(d){
  AL <- matrix(NA, d, d)
  ind <- t(which(upper.tri(AL), arr.ind=TRUE))
  qind=length(ind[1,])
  return(list(ind = ind, qind = qind))
}

# for reordering a parameter vector created from a matrix using utils::combn
#ltheta is the length of the parameter vector defining the matrix
combparam2uppertriorder_matrix <- function(ltheta){
  # l(theta) = p-1 + (p-1) * (p-2)/2
  # = p-1 + (p^2 -p -2p + 2)/2
  # 2l = 2p-2 + p^2 - 3p + 2
  # 0 = p^2  - p - 2l
  # p = (1 pm sqrt(1 - 4*(-2l)))/2
  #   = (1 pm sqrt(8l + 1)) /2
  #   = (1+sqrt(8l + 1) / 2
  p <- (1+sqrt(8 * ltheta + 1)) / 2
  sp <- p-1

  #upper.tri ordering
  AL <- matrix(NA, sp, sp)
  diag(AL) <- 1:sp
  AL[upper.tri(AL)] <- sp + (1:sum(upper.tri(AL)))
  AL[lower.tri(AL)] <- t(AL)[lower.tri(AL)]

  #comb ordering
  inthetaorder <- c(diag(AL),
                    AL[t(combn(1:(p-1), 2))]) #each column is the own dimension, each row single element to extract

  stopifnot(length(inthetaorder) == ltheta) #check that length produced is correct
  stopifnot(all(inthetaorder[inthetaorder] == seq(1, length.out = ltheta))) #check that ordering index is correct
  return(inthetaorder)
}

