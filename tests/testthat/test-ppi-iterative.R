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

  hardcodedestimate <- estimator1(model$sample, acut, incb = TRUE, beta = model$beta)
  
  expect_absdiff_lte_v(out$est$paramvec, hardcodedestimate$est$paramvec, 0.0001 * out$SE$paramvec) #proxy for optimisation flatness
  expect_absdiff_lte_v(out$est$paramvec, model$theta, 2 * out$SE$paramvec)
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
  expect_equal(est_cppad$est$paramvec, estimator$est$paramvec, tolerance = 1E-3, ignore_attr = "names")
})

test_that("Iterative solver works on microbiome data with outliers, alr", {
  skip("next calculation, the cppad estimate, takes hours")
  list2env(ppi_microbiomedata_TCAP(), globalenv())
  
  #hardcoded estimate:
  est_hardcoded=ppi(Y = propreal,
         method = "hardcoded", trans = "alr",
         paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0))

system.time({est_iterative=ppi(Y = propreal,
         method = "iterative", trans = "alr",
         paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
         bdrythreshold = 1E-20, shiftsize = 1E-20,
         control = list(maxit = 1E5, tol = 1E-20 * 100))})
# the defaults for Rcgmin are meaning the estimate takes too long to converge!
# from the default starting parameter vec it takes many more iterations than normal to get to converge to the correct result.
# With current tolerance of 1E-20*100, it took 4864 seconds (1.3+ hours) and 800000 grad evaluations, and still hadn't reached the tolerance.
# however it close enough to to the hardcoded estimate
expect_equal(est_hardcoded$est$paramvec, est_iterative$est$paramvec, tolerance = 1E-3)
})

test_that("cppad_search gives similar result to cppad_closed", {
  set.seed(354)
  m <- ppi_egmodel(100, maxden = 4)
  tapes <- buildsmotape("sim","sqrt", "sph", "ppi",
                        ytape = rep(1/m$p, m$p),
                        usertheta = rep(NA, length(m$theta)),
                        divweight = "minsq", acut = 0.1)
  estsearch <- cppad_search(tapes$smotape, m$theta *0 + 1, m$sample, control = list(tol = 1E-14, maxit = 2000))
  estclosed <- cppad_closed(tapes$smotape, m$sample)
  expect_equal(estsearch$est, estclosed$est, tolerance = 1E-3, ignore_attr = "names")
  expect_equal(estsearch$SE, estclosed$SE, tolerance = 1E-3)
})

test_that("cppad_search with weights gives similar result to cppad_closed", {
  set.seed(124)
  m <- ppi_egmodel(100, maxden = 4)
  tapes <- buildsmotape("sim","sqrt", "sph", "ppi",
                        ytape = rep(1/m$p, m$p),
                        usertheta = rep(NA, length(m$theta)),
                        divweight = "minsq", acut = 0.1)
  w <- runif(100)
  estsearch <- cppad_search(tapes$smotape, m$theta *0 + 1, m$sample, control = list(tol = 1E-12, maxit = 1000), w = w)
  estclosed <- cppad_closed(tapes$smotape, m$sample, w = w)
  expect_equal(estsearch$est, estclosed$est, tolerance = 1E-2, ignore_attr = "names")
  expect_equal(estsearch$SE, estclosed$SE)
})

test_that("cppad_search output value matches smvalues_tape result", {
  set.seed(1234)
  m <- ppi_egmodel(1000, maxden = 4)
  tapes <- buildsmotape("sim","sqrt", "sph", "ppi",
                        ytape = rep(1/m$p, m$p),
                        usertheta = rep(NA, length(m$theta)),
                        divweight = "minsq", acut = 0.1)

  est <- cppad_search(tapes$smotape, m$theta *0 + 1, m$sample, control = list(tol = 1E-12, maxit = 1000))

  smvalues <- smvalues_tape_wsum(tapes$smotape, m$sample, est$est)
  expect_equal(est$value, smvalues$obj/nrow(m$sample))
  expect_equal(est$sqgradsize, sum(smvalues$grad^2))
})

