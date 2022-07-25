test_that("ppi Dirichlet with cW gives correct params on simulated data, with two outliers. p=3", {
  set.seed(1273)
  m <- ppi_egmodel(1000, maxden = 4)
  outlier1 <- c(0.9, 0.9, 0.01)
  outlier1 <- outlier1/sum(outlier1)
  outlier2 <- c(0.9, 0.1, 0.01)
  outlier2 <- outlier2/sum(outlier2)
  m$sample <- rbind(m$sample, outlier1, outlier2)

  #non-robust estimates
  est_norobust <- ppi(m$sample, AL = 0, bL = 0, acut=0.1, method = "direct", trans = "sqrt", bdryweight = "minsq")
  est_norobust2 <- ppi(m$sample, AL = 0, bL = 0, acut=0.1, method = "direct", trans = "sqrt", bdryweight = "minsq", cW = 0)

  expect_equal(est_norobust2$theta, est_norobust$theta)

  #robust
  est_robust1 <- ppi(m$sample, AL = 0, bL = 0, acut=0.1, method = "direct", trans = "sqrt", bdryweight = "minsq", cW = 1E-1)

  rmse <- function(v1, v2){sqrt(mean((v1 - v2)^2))}
  expect_lt(rmse(m$theta, est_norobust$theta),
            rmse(m$theta, est_robust1$theta))
})

test_that("Robustness runs for direct and cppad methods", {
  set.seed(1273)
  m <- ppi_egmodel(50, maxden = 4)  # at 20 got singularities
  outlier1 <- c(0.9, 0.9, 0.01)
  outlier1 <- outlier1/sum(outlier1)
  outlier2 <- c(0.9, 0.1, 0.01)
  outlier2 <- outlier2/sum(outlier2)
  m$sample <- rbind(m$sample, outlier1, outlier2)

  #Ralr
  out <- ppi(m$sample, bL = 0, betap = -0.5, method = "direct", trans = "alr", cW = 1E-1)
  expect_gt(out$optim$fpevals, 1)

  #dir minsq
  out <- ppi(m$sample, AL = 0, bL = 0, method = "direct",
             trans = "sqrt", acut = 0.1, bdryweight = "minsq", cW = 1E-1)
  expect_gt(out$optim$fpevals, 1) #errors
  #dir prodsq
  out <- ppi(m$sample, AL = 0, bL = 0, method = "direct",
             trans = "sqrt", acut = 0.1, bdryweight = "prodsq", cW = 1E-1)
  expect_gt(out$optim$fpevals, 1) #errors

  # estimator1 bL=0
  out <- ppi(m$sample, bL = 0, beta = m$beta, acut=0.1, method = "direct", trans = "sqrt", bdryweight = "minsq", cW = 1E-2)
  expect_gt(out$optim$fpevals, 1)
  # estimator1 bL=0 prodsq
  out <- ppi(m$sample, bL = 0, beta = m$beta,
             acut=0.1, method = "direct", trans = "sqrt", bdryweight = "prodsq",
             cW = 1E-2)
  expect_gt(out$optim$fpevals, 1)

  # estimator1 bL!=0
  out <- ppi(m$sample, beta = m$beta, acut=0.1, method = "direct", trans = "sqrt", bdryweight = "minsq", cW = 1E-2)
  expect_gt(out$optim$fpevals, 1)
  # estimator1 bL!=0 prodsq
  out <- ppi(m$sample, beta = m$beta,
             acut=0.1, method = "direct", trans = "sqrt", bdryweight = "prodsq", cW = 1E-2)
  expect_gt(out$optim$fpevals, 1)
  # estimatorall1 betap fixed
  out <- ppi(m$sample, betap = tail(m$beta, 1),
             acut=0.1, method = "direct", trans = "sqrt", bdryweight = "minsq",
             cW = 1E-2)
  expect_gt(out$optim$fpevals, 1)

  # estimatorall1 betap fitted
  out <- ppi(m$sample,
             acut=0.1, method = "direct", trans = "sqrt", bdryweight = "minsq",
             cW = 1E-2)
  expect_gt(out$optim$fpevals, 1)

  # cppad
  out <- ppi(m$sample,
             acut=0.1, method = "cppad", trans = "sqrt", bdryweight = "minsq",
             cW = 1E-2)
  expect_gt(out$optim$fpevals, 1)
})
