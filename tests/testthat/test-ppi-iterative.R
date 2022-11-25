# test of iterative solver for ppi (compare to hardcoded):
## simple
## weighted fit
## microbiomfit without outliers
## microbiomfit with outliers

test_that("ppi iterative solve match estimator1 minsq with fixed beta for ppi_egmodel", {
  set.seed(123)
  model <- ppi_egmodel(1000, maxden = 4)

  acut = 0.1
  out <- ppi(model$sample, paramvec = ppi_paramvec(beta = model$beta),
            method = "iterative",
            bdrythreshold = 0,
            trans = "sqrt", divweight = "minsq", acut = acut)

  directestimate <- estimator1(model$sample, acut, incb = TRUE, beta = model$beta0)
  
  expect_absdiff_lte_v(out$est$paramvec, directestimate$est$paramvec, 0.0001 * out$SE$paramvec) #proxy for optimisation flatness
  expect_absdiff_lte_v(out$est$paramvec, model$theta, 2 * out$SE$paramvec)
})


