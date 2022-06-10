test_that("to then from PPI param vector returns the same params", {
  set.seed(124)
  p = 5
  ALs <- rsymmetricmatrix(p-1, 0, 10)
  bL <- runif(p-1)
  beta <- runif(p, -1, 0)
  paramvec <- toPPIparamvec(ALs,  bL, beta)
  params2 <- fromPPIparamvec(paramvec, p)
  expect_equal(params2$ALs, ALs)
  expect_equal(params2$bL, bL)
  expect_equal(params2$beta, beta)
})

test_that("from PPI param vector order", {
  set.seed(124)
  p = 5
  bL <- runif(p-1)
  beta <- runif(p, -1, 0)
  paramats <- fromPPIparamvec(c(seq(1:((p-1) + (p-1) * (p-2)/2)), bL, beta), p)

  #alternative
  sp <- p-1
  x=c(1:sp)
  ind=combn(x, 2, FUN = NULL, simplify = TRUE)
  qind=length(ind[1,])
  mat=matrix(NA,sp,sp)
  vals <- seq(1:((p-1) + (p-1) * (p-2)/2))
  for (j in 1:sp)
  {
    mat[j,j]=vals[j]
  }
  for (j in 1:qind)
  {
    mat[ind[1,j],ind[2,j]]=vals[sum(j,sp)]
    mat[ind[2,j],ind[1,j]]=vals[sum(j,sp)]
  }

  expect_false(all(mat == paramats$ALs))

})
