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
  ALs <- rsymmetricmatrix(p-1, -4, 4)
  bL <- matrix(0, nrow = p-1)
  beta <- c(-0.7, -0.8, -0.3, 0, 0)
  set.seed(13456) #this seed generates samples that are representative-enough for estimatorlog_ratio() to give close estimates
  sim <- rhybrid(1000, p, beta, ALs, bL, 4)
  prop <- sim$samp3
  # prop %>% as_tibble() %>% tidyr::pivot_longer(everything()) %>% ggplot() + facet_wrap(vars(name)) + geom_freqpoly(aes(x=value))

  #calculate robust estimates
  cW=0.1
  est1=windham_diff(prop,cW,ALs,bL,beta, ind_weightA = c(0,0,0,1), originalcorrectionmethod = TRUE)
  expect_equal(est1$est$ALs, ALs, tolerance = 1E0)
  expect_equal(est1$est$beta, beta, tolerance = 1E-1)

  est2=windham_diff(prop,cW,ALs,bL,beta, ind_weightA = c(0,0,0,1), originalcorrectionmethod = FALSE)
  expect_equal(est1$est$theta, est2$est$theta)
})
