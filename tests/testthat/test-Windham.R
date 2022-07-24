test_that("windam_raw gives correct params on simulated data, with two outliers. p=3", {
  skip("tested by all robust estimator tests")
  set.seed(1273)
  m <- ppi_egmodel(1000, maxden = 4)
  outlier1 <- c(0.9, 0.9, 0.01)
  outlier1 <- outlier1/sum(outlier1)
  outlier2 <- c(0.9, 0.1, 0.01)
  outlier2 <- outlier2/sum(outlier2)
  m$sample <- rbind(m$sample, outlier1, outlier2)

  #check non-robust estimates, excluding the outliers
  est_simple <- estimatorall1(m$sample[1:1000, ], acut=0.1)
  #get non-robust estimates with the outliers
  est_simple_outlier <- estimatorall1(m$sample, acut=0.1)

  #calculate robust estimates
  inWW <- ppi_cppad_thetaprocessor(m$p, AL = TRUE, bL = FALSE, beta = FALSE) #all dimensions have negative beta.
  ppildenfun <- function(sample, theta){
    ppiparmats <- fromPPIparamvec(theta)
    logden <- dppi(sample, ppiparmats$beta, ppiparmats$ALs, ppiparmats$bL)
    return(logden)
  }
  ppiestimator <- function(Y, starttheta, isfixed, w){
    estimatorall1(Y, acut = 0.1, w = w)$estimator1
  }

  isfixed <- ppi_cppad_thetaprocessor(m$p, AL=FALSE, bL = FALSE, betaL = FALSE, betap = FALSE)
  est <- windham_raw(prop = m$sample,
                     cW = 0.1 * inWW,
                     ldenfun = ppildenfun,
                     estimatorfun = ppiestimator,
                     starttheta = m$theta * 0,
                     isfixed = isfixed,
                     originalcorrectionmethod = TRUE)

  expect_equal(est$theta, est_simple$estimator1, tolerance = 0.1, ignore_attr = TRUE)
  #below checks that the non-robust estimate with outliers is much different to the robust estimate
  expect_error(expect_equal(est$theta, est_simple_outlier$estimator1, tolerance = 0.1, ignore_attr = TRUE))
})

test_that("windam_diff with estimatorall1 gives correct params on simulated, no outlier, data. p=3", {
  set.seed(1273)
  p = 3
  ALs <- exp(rsymmetricmatrix(p-1, -4, 4))
  bL <- matrix(0, nrow = p-1)
  beta <- c(-0.8, -0.3, 0)
  set.seed(1345) #this seed generates samples that are representative-enough for estimatorlog_ratio() to give close estimates
  prop <- rppi(1000, p, beta, ALs, bL, 4)$samp3

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
  sim <- rppi(1000, p, beta, ALs, bL, 4)
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
