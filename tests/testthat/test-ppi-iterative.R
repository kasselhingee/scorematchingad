# test of iterative solver for ppi (compare to hardcoded):
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


test_that("Simulated weights for ppi with minsq match itself and estimatorall1", {
  set.seed(1243)
  m <- ppi_egmodel(1000, maxden = 4)
  #simulate weights
  set.seed(134)
  vw <- virtualweights(m$sample)
  acut = 0.01

  out_sim <- ppi(vw$newY, 
            method = "iterative",
            bdrythreshold = 0,
            trans = "sqrt", divweight = "minsq", acut = 0.01,
            control = list(tol = 1E-20, maxit = 2000))
  out_dir <- ppi(m$sample, 
            method = "iterative",
            bdrythreshold = 0,
            trans = "sqrt", divweight = "minsq", acut = 0.01,
            w = vw$w,
            control = list(tol = 1E-20, maxit = 2000))

  expect_equal(out_sim[!(names(out_sim) %in% c("counts", "SE"))], 
     out_dir[!(names(out_sim) %in% c("counts", "SE"))],
     tolerance = 1E-3)

  directestimate <- estimatorall1(m$sample, acut = 0.01, w = vw$w)
  expect_equal(directestimate$estimator1, out_sim$est$paramvec, ignore_attr = TRUE)
  expect_equal(directestimate$estimator1, out_dir$est$paramvec, ignore_attr = TRUE)
})

test_that("Iterative solver works on microbiome data without outliers", {
  list2env(ppi_microbiomedata_prep1(), globalenv())
  
  #hardcoded estimate:
  estimator= estimator1(propreal, 0.01,1, beta0, computeSE = TRUE)
  #check it matches iterative ppi()
  est_cppad <- ppi(Y = propreal, acut = 0.01,
                   method = "closed",
                   trans = "sqrt", divweight = "minsq",
                   bdrythreshold = 1E-15, shiftsize = 1E-10,
                   paramvec = ppi_paramvec(beta = beta0),
                   control = list(maxit = 100000, tol = 1E-20))
  expect_equal(est_cppad$est$paramvec, estimator$est$paramvec, tolerance = 1E-3)

 
})
