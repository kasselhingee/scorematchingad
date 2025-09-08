test_that("estimating matrix von Mises-Fisher get close to correct", {
  set.seed(1)
  M <- matrix(rnorm(4*2), nrow = 4, ncol = 2)
  set.seed(3)
  samp <- replicate(100, rstiefel::rmf.matrix(M))
  
  set.seed(4)
  Mother <- matrix(rnorm(nrow(M)*ncol(M)),
                   nrow = nrow(M), ncol = ncol(M))
  tape_smd("Steifel", nrow = nrow(M), ncol = ncol(M),
           ll = "Stiefel_MF",
           ytape = vec(samp[,,1]),
           usertheta = vec(Mother))
})