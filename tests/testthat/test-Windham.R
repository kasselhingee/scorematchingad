test_that("windam_diff gives correct params on simulated, no outlier, data. p=3", {
  set.seed(1273)
  p = 3
  ALs <- exp(rsymmetricmatrix(p-1, -4, 4))
  bL <- matrix(0, nrow = p-1)
  beta <- c(-0.8, -0.3, 0)
  set.seed(1345) #this seed generates samples that are representative-enough for estimatorlog_ratio() to give close estimates
  prop <- rhybrid(1000, p, beta, ALs, bL, 4)$samp3

  #check non-robust estimates
  est_unload <- estimatorlog_weight(prop, betap = beta[p], weightW = rep(1, nrow(prop)))
  # fromPPIparamvec(est_unload$ppi, p)$ALs #fairly terrible at the AL
  # fromPPIparamvec(est_unload$ppi, p)$beta #pretty good at beta

  #calculate robust estimates
  cW=0.1
  est1=windham_diff(prop,cW,ALs,bL,beta, ind_weightA = c(0,0), originalcorrectionmethod = TRUE)
  expect_equal(est1$est$ALs, ALs, tolerance = 1)
  expect_equal(est1$est$beta, beta, tolerance = 1E-1)

  est2=windham_diff(prop,cW,ALs,bL,beta, ind_weightA = c(0,0), originalcorrectionmethod = FALSE)
  expect_equal(est2$est$ALs, ALs, tolerance = 1)
  expect_equal(est2$est$beta, beta, tolerance = 1E-1)

  expect_equal(est2$est$theta, est1$est$theta)
})

test_that("windam_diff gives correct params on simulated, no outlier, data. p = 5", {
  set.seed(1273)
  p = 5
  ALs <- exp(rsymmetricmatrix(p-1, -4, 4))
  bL <- matrix(0, nrow = p-1)
  beta <- c(-0.7, -0.8, -0.3, 0, 0)
  set.seed(13456) #this seed generates samples that are representative-enough for estimatorlog_ratio() to give close estimates
  prop <- rhybrid(1000, p, beta, ALs, bL, 35)$samp3

  #check non-robust estimates
  est_unload <- estimatorlog_weight(prop, betap = beta[p], weightW = rep(1, nrow(prop)))
  fromPPIparamvec(est_unload$ppi, p)$ALs #fairly terrible at the AL
  fromPPIparamvec(est_unload$ppi, p)$beta #pretty bad at beta too: -10 and -15 isn't even defined!?

  est_unload_cppad <- ppi_cppad(prop, bL = 0, betap=beta[p], man = "Ralr", weightname = "ones")
  est_unload_cppad$est$beta #pretty bad at beta still: -4 and -2 isn't even defined!?

  #calculate robust estimates
  cW=0.1
  est1=windham_diff(prop,cW,ALs,bL,beta, ind_weightA = c(0,0,0,1), originalcorrectionmethod = TRUE)
  expect_equal(est1$est$ALs, ALs, tolerance = 1E-1)

  est2=windham_diff(prop,cW,ALs,bL,beta, ind_weightA = c(0,0,0,1), originalcorrectionmethod = FALSE)
  expect_equal(est2$est$ALs, ALs, tolerance = 1E-1)

  est_unload <- estimatorlog_weight(prop, betap = beta[p], weightW = rep(1, nrow(prop)))
  fromPPIparamvec(est_unload$ppi, p)$ALs
  ALs

  est_cppad <- ppi_cppad(prop, bL = rep(0, p-1), betap = beta[p], man = "Ralr", weightname = "ones",
                    control = list(tol = 1E-10))
})
