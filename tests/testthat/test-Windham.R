test_that("windam_raw gives correct params on simulated data, with two outliers. p=3", {
  skip_on_cran("only extra checks are the variable cW ones")
  set.seed(1273)
  m <- ppi_egmodel(1000, maxden = 4)
  outlier1 <- c(0.9, 0.9, 0.01)
  outlier1 <- outlier1/sum(outlier1)
  outlier2 <- c(0.9, 0.1, 0.01)
  outlier2 <- outlier2/sum(outlier2)
  m$sample <- rbind(m$sample, outlier1, outlier2)

  #check non-robust estimates, excluding the outliers
  est_simple <- ppi(m$sample[1:1000, ], acut=0.1, method = "direct", trans = "sqrt", bdryweight = "minsq")
  #get non-robust estimates with the outliers
  est_simple_outlier <- ppi(m$sample, acut=0.1, method = "direct", trans = "sqrt", bdryweight = "minsq")

  #calculate robust estimates
  ppildenfun <- function(sample, theta){
    ppiparmats <- cdabyppi:::fromPPIparamvec(theta)
    logden <- dppi(sample, ppiparmats$beta, ppiparmats$ALs, ppiparmats$bL)
    return(logden)
  }
  ppiestimator <- function(Y, starttheta, isfixed, w){
    ppi(Y, acut = 0.1, w = w, method = "direct", trans = "sqrt", bdryweight = "minsq")$estimator1
  }

  isfixed <- cdabyppi:::ppi_cppad_thetaprocessor(m$p, AL=FALSE, bL = FALSE, betaL = FALSE, betap = FALSE)
  est <- cdabyppi:::windham_raw(prop = m$sample,
                     cW = 0.1 * cdabyppi:::ppi_cppad_thetaprocessor(m$p, AL = TRUE, bL = FALSE, beta = FALSE), #all dimensions have negative beta.
                     ldenfun = ppildenfun,
                     estimatorfun = ppiestimator,
                     starttheta = m$theta * 0,
                     isfixed = isfixed,
                     originalcorrectionmethod = TRUE)

  # variable c, expect estimates to be different
  cW <- cdabyppi:::ppi_cppad_thetaprocessor(m$p, AL = matrix(c(0.1, 1E-3, 1E-3, 0.1), nrow = 2, ncol = 2),
                                 bL = 0, beta = 0)
  expect_error(cdabyppi:::windham_raw(prop = m$sample,
                     cW = cW,
                     ldenfun = ppildenfun,
                     estimatorfun = ppiestimator,
                     starttheta = m$theta * 0,
                     isfixed = isfixed,
                     originalcorrectionmethod = TRUE)) #error because original correction method doesn't cope with variable non-zero cW elements
  est_varcW <-  cdabyppi:::windham_raw(prop = m$sample,
                     cW = cW,
                     ldenfun = ppildenfun,
                     estimatorfun = ppiestimator,
                     starttheta = m$theta * 0,
                     isfixed = isfixed,
                     originalcorrectionmethod = FALSE)

  expect_equal(est$theta, est_simple$estimator1, tolerance = 0.1, ignore_attr = TRUE)
  #below checks that the non-robust estimate with outliers is much different to the robust estimate
  expect_error(expect_equal(est$theta, est_simple_outlier$estimator1, tolerance = 0.1, ignore_attr = TRUE))

  # expect that the different cW values would lead to different estimates
  expect_gt(mean(abs(est$theta - est_varcW$theta)), 10)
})

test_that("windam_diff with gives correct params on simulated, no outlier, data. p=3", {
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

  rmse <- function(v1, v2){sqrt(mean((v1 - v2)^2))}
  expect_gt(rmse(beta, est_unload$ppi[6:8]), rmse(beta, est1$est$beta))
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


