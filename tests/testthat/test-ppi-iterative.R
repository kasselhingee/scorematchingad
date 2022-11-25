# test of iterative solver for ppi (compare to hardcoded):
## microbiomfit with outliers

test_that("ppi iterative solve match estimator1 minsq with fixed beta for ppi_egmodel", {
  set.seed(123)
  model <- ppi_egmodel(1000, maxden = 4)

  acut = 0.1
  out <- ppi(model$sample, paramvec = ppi_paramvec(beta = model$beta),
            method = "iterative",
            bdrythreshold = 0,
            trans = "sqrt", divweight = "minsq", acut = acut)

  hardcodedestimate <- estimator1(model$sample, acut, incb = TRUE, beta = model$beta0)
  
  expect_absdiff_lte_v(out$est$paramvec, hardcodedestimate$est$paramvec, 0.0001 * out$SE$paramvec) #proxy for optimisation flatness
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
            control = list(tol = 1E-15, maxit = 2000))
  out_dir <- ppi(m$sample, 
            method = "iterative",
            bdrythreshold = 0,
            trans = "sqrt", divweight = "minsq", acut = 0.01,
            w = vw$w,
            control = list(tol = 1E-15, maxit = 2000))

  expect_equal(out_sim[!(names(out_sim) %in% c("info", "SE"))], 
     out_dir[!(names(out_sim) %in% c("info", "SE"))],
     tolerance = 1E-2)

  hardcodedestimate <- estimatorall1(m$sample, acut = 0.01, w = vw$w)
  expect_equal(hardcodedestimate$estimator1, out_sim$est$paramvec, tolerance = 1E-2, ignore_attr = TRUE)
})

test_that("Iterative solver works on microbiome data without outliers", {
  list2env(ppi_microbiomedata_cleaned_TCAP(), globalenv())
  
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

test_that("Iterative solver works on microbiome data with outliers, alr", {
  #skip("next calculation, the cppad estimate, takes hours")
  list2env(ppi_microbiomedata_TCAP(), globalenv())
  
  #hardcoded estimate:
  est_hardcoded=ppi(Y = propreal,
         method = "hardcoded", trans = "alr",
         paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0))

system.time({est_cppad=ppi(Y = propreal,
         method = "closed", trans = "alr",
         paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
         bdrythreshold = 1E-20, shiftsize = 1E-20,
         control = list(maxit = 1E5, tol = 1E-20 * 100))})
# the defaults for Rcgmin are meaning the estimate takes too long to converge!
# from the default starting parameter vec it takes many more iterations than normal to get to converge to the correct result. In this case of default tolerance, then 43961 evaluations of the gradient.
# With current tolerance of 1E-20*n, then 4559 seconds (1.3 hours), 81479 grad evaluations.
expect_equal(est_hardcoded$est$paramvec, est_cppad$est$paramvec, tolerance = 1E-3)
})
