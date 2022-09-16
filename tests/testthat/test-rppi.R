# test that rppi old matches old rppi new

#sample size
n=5000

#parameters for the PPI model
m <- ppi_egmodel(2)

test_that("current PPI simulation method gives samples with similar empirical density estimates as the original simulation method", {
  skip_on_cran() #accuracy of the method is tested by all the estimators
  # simulate using old method
  time_historic <- system.time(samp2 <- cdabyppi:::rppi_singly(n,m$p,m$beta0,m$ALs,m$bL,4))
  H <- ks::Hpi(samp2$samp3[, -m$p])
  kde_historic <- ks::kde(samp2$samp3[, -m$p], H)
  #simulate sample from PPI model
  time_current <- system.time(samp1 <- cdabyppi:::rppi(n,m$p,m$beta0,m$ALs,m$bL,4))
  H <- ks::Hpi(samp1$samp3[, -m$p])
  kde_current <- ks::kde(samp1$samp3[, -m$p], H)

  testresult <- ks::kde.test(samp2$samp3[, -m$p], samp1$samp3[, -m$p])
  expect_gte(testresult$pvalue, 0.01)
  time_ratio <- time_current/time_historic
  expect_lt(time_ratio[["user.self"]], 0.05)
})

test_that("rppi() is fixed by set.seed()", {
  m <- ppi_egmodel(2)

  set.seed(3212)
  Y1 <- rppi(100, m$p, m$beta, m$AL, m$bL, 4)

  set.seed(3212)
  Y2 <- rppi(100, m$p, m$beta, m$AL, m$bL, 4)

  expect_equal(Y1, Y2)
})
