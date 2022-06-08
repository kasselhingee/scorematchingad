# test that rhybrid old matches old rhybrid new

#sample size
n=5000

#parameters for the PPI model
m <- sec2_3model(2)

test_that("current PPI simulation method gives samples with similar empirical density estimates as the original simulation method", {
  # simulate using old method
  time_historic <- system.time(samp2 <- cdabyppi:::rhybrid_singly(n,m$p,m$beta0,m$ALs,m$bL,4))
  H <- ks::Hpi(samp2$samp3[, -p])
  kde_historic <- ks::kde(samp2$samp3[, -p], H)
  #simulate sample from PPI model
  time_current <- system.time(samp1 <- cdabyppi:::rhybrid(n,m$p,m$beta0,m$ALs,m$bL,4))
  H <- ks::Hpi(samp1$samp3[, -p])
  kde_current <- ks::kde(samp1$samp3[, -p], H)

  testresult <- ks::kde.test(samp2$samp3[, -p], samp1$samp3[, -p])
  expect_gte(testresult$pvalue, 0.01)
  time_ratio <- time_current/time_historic
  expect_lt(time_ratio[["user.self"]], 0.02)
})

