test_that("PPI ALR hardcoded estimate has low smgrad values", {
  mnongamma <- ppi_egmodel(1)
  theta <- ppi_paramvec(beta = c(-0.95, -0.9, 0.5), AL = mnongamma$AL, bL = 0)
  set.seed(1234)
  Ycts <- rppi(1000, paramvec = theta)
  dsample <- round(Ycts * 100)/ 100
  dsample[, 3] <- 1 - rowSums(dsample[, 1:2])
  colMeans(dsample == 0)
  mean(apply(dsample, 1, min) == 0)  #0.96

  usertheta = ppi_paramvec(p = ncol(dsample),
                           bL = 0,
                           betap = tail(theta, 1))
  est_hardcoded <- ppi(dsample, paramvec = usertheta, trans = "alr", method = "hardcoded")

  ppi_smvalues(dsample, paramvec = usertheta, evalparam = est_hardcoded$est$paramvec, trans = "alr")

  tapes <- buildsmotape("Ralr", "ppi",
                utape = rep(1/mnongamma$p, mnongamma$p),
                usertheta = usertheta)

  hardcodedvals <- tape_smvalues(tapes$smotape, 
         xmat = dsample,
         pmat = t_ut2f(usertheta, est_hardcoded$est$paramvec))
  modelvals <- ppi_cppad_values(dsample,
         stheta = theta,
         isfixed = t_u2i(usertheta),
         man = "Ralr",
         hsqfun = "ones", 
         acut = 1)

  expect_lt(hardcodedvals$obj, modelvals$obj) #because for any given sample the estimate would be better than the true value

  expect_lt_v(abs(hardcodedvals$grad), abs(modelvals$grad)) #because the estimate will be better for any given sample
  expect_lt_v(abs(hardcodedvals$grad), rep(1E-10, length(hardcodedvals$grad)))
})


