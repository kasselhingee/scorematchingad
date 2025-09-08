test_that("estimating matrix von Mises-Fisher get close to correct", {
  M <- matrix(rnorm(4*2), nrow = 4, ncol = 2)
  set.seed(3)
  samp <- replicate(100, rstiefel::rmf.matrix(M))
})