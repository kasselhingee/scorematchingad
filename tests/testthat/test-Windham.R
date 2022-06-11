test_that("windam_diff gives correct params on simulated, no outlier, data", {
  set.seed(1273)
  p = 5
  ALs <- exp(rsymmetricmatrix(p-1, -4, 4))
  bL <- matrix(0, nrow = p-1)
  beta <- c(-0.7, -0.8, -0.3, 0, 0)
  prop <- rhybrid(1000, p, beta, ALs, bL, 35)$samp3

  #calculate robust estimates
  cW=0.1
  est1=windham_diff(prop,cW,ALs,bL,beta, ind_weightA = c(0,0,0,1))
  expect_equal(est1$ALs_est, ALs, tolerance = 1E-1)

  est_unload <- estimatorlog_weight(prop, betap = beta[p], weightW = rep(1, nrow(prop)))
  fromPPIparamvec(est_unload$ppi, p)$ALs
  est1$ALs_est

  est_cppad <- ppi_cppad(prop, bL = rep(0, 3-1), betap = m$beta0[p], man = "Ralr", weightname = "ones",
                    control = list(tol = 1E-10))
})
