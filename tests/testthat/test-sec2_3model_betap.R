####model in Section 2.3 with beta=(-0.8,-0.8,?) here a_c is set to 20%####
# edited from 'modelA.R' in the original code

#### Setup ####
list2env(ppi_egmodel(1), globalenv())
theta <- c(diag(ALs), ALs[upper.tri(ALs)], bL)



#### Tests ####
test_that("Score1ac estimator estimates beta0[0] and other consistently with cppad version", {
  #sample size
  n=100
  beta0=matrix(-0.8,p,1)
  beta0[p]=-0.5

  #simulate sample from PPI model
  set.seed(321)
  samp1=cdabyppi:::rppi(n,p,beta0,ALs,bL,4)
  samp3=samp1$samp3
  theta <- cdabyppi:::ppi_paramvec(p,
              AL = ALs, bL = drop(bL), beta = drop(beta0))

  ####Score1ac estimator##

  #a_c for h function:
  acut=1.6e-02

  #calculate scoring estimate for full model
  estimator=cdabyppi:::estimatorall1(samp3,acut, betap = NULL)
  estimate1all=estimator$estimator1

  # Get SE from CppAD methods
  intheta <- ppi_paramvec(p)
  tapes <- buildsmotape("sphere", "ppi",
                        samp3[1, ], intheta,
                        weightname = "minsq",
                        acut = acut)
  SE <- cppadSE(tapes$smotape, estimate1all, samp3)

  expect_absdiff_lte_v(estimate1all, theta, 3 * SE)

  # compare to ppi via cppad
  expect_lt(sum(smobjgrad(tapes$smotape, estimate1all, samp3)^2), 1E-14)
  est2 <- ppi(samp3, trans = "sqrt", bdryweight = "minsq", acut = acut, bdrythreshold = 1E-20, control = list(tol = 1E-15), method = "cppad", w = NULL)
  expect_equal(est2$est$paramvec, drop(estimate1all), tolerance = 1E-2) #within 1% of each other roughly
  cdabyppi:::expect_absdiff_lte_v(drop(estimate1all), est2$est$paramvec, 2*est2$SE$paramvec)
})

test_that("Score1ac estimator can estimate beta0[1:(p-1)] for beta0[p] larger than -0.5", {
  #sample size
  n=100
  beta0=matrix(-0.8,p,1)
  beta0[p] = 5
  theta <- cdabyppi:::ppi_paramvec(p,
                                               AL = ALs, bL = drop(bL), beta = drop(beta0))

  #simulate sample from PPI model
  set.seed(124)
  samp1=cdabyppi:::rppi(n,p,beta0,ALs,bL,4)
  samp3=samp1$samp3

  stopifnot(samp1$maxden <= 4)
  maxden <- 4

  ####Score1ac estimator##

  #a_c for h function:
  acut=1.6e-02

  #calculate scoring estimate for full model, beta[p] correctly fixed
  estimator=cdabyppi:::estimatorall1(samp3,acut, betap = beta0[p])
  estimate1all=estimator$estimator1

  # SE from cppad
  intheta <- cdabyppi:::ppi_paramvec(p, betap = beta0[p])
  tapes <- buildsmotape("sphere", "ppi",
                        samp3[1, ], intheta,
                        weightname = "minsq",
                        acut = acut)
  SE <- cppadSE(tapes$smotape, estimate1all, samp3)

  #3*SE bounds
  expect_absdiff_lte_v(estimate1all, theta[is.na(intheta)], 3 * SE)
})

test_that("Score1ac estimator can estimate beta0[1:(p-1)] for beta0[p] large but misspecified", {
  #sample size
  n=100
  beta0=matrix(-0.8,p,1)
  beta0[p]= 5
  theta <- cdabyppi:::ppi_paramvec(p,
                                               AL = ALs, bL = drop(bL), beta = drop(beta0))

  #simulate sample from PPI model
  samp1=cdabyppi:::rppi(n,p,beta0,ALs,bL,4)
  samp3=samp1$samp3

  stopifnot(samp1$maxden <= 4)
  maxden <- 4

  ####Score1ac estimator##

  #a_c for h function:
  acut=1.6e-02

  #calculate scoring estimate for full model, beta[p] correctly fixed
  estimator=cdabyppi:::estimatorall1(samp3,acut, betap = -0.5)
  estimate1all=estimator$estimator1

  # SE from cppad
  intheta <- cdabyppi:::ppi_paramvec(p, betap = -0.5)
  tapes <- buildsmotape("sphere", "ppi",
                        samp3[1, ], intheta,
                        weightname = "minsq",
                        acut = acut)
  SE <- cppadSE(tapes$smotape, estimate1all, samp3)

  #3*SE bounds
  expect_absdiff_lte_v(estimate1all, theta[is.na(intheta)], 3 * SE)
})
