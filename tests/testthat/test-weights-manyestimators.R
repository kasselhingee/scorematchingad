set.seed(1234)
m <- ppi_egmodel(1000, maxden = 4)
#simulate weights
set.seed(134)
vw <- virtualweights(m$sample)
acut = 0.1

test_that("estimator2 and SE matches for simulated weights", {
  est_sim <- estimator2(vw$newY, acut = acut, incb = 1, beta0 = m$beta0)
  est_direct <- estimator2(m$sample, acut = acut, incb = 1, beta0 = m$beta0, w = vw$w)
  expect_equal(est_direct, est_sim)

  estSE_sim <- estimator2SE(vw$newY,acut,est_sim$estimator2, est_sim$W_est, incb = 1, beta0 = m$beta0)
  estSE_direct <- estimator2SE(m$sample, acut = acut, est_sim$estimator2, est_sim$W_est, incb = 1, beta0 = m$beta0, w = vw$w)
  expect_equal(estSE_direct, estSE_sim)
})

test_that("ppi_mmmm matches for simulated weights", {
  #simulate sample from the multinomial PPI model:
  x=matrix(0,nrow(m$sample),m$p)
  for (j in 1:nrow(m$sample))
  {
    x[j,]=rmultinom(1,2000,prob=m$sample[j,])
  }
  vw <- virtualweights(x)
  est_sim <- ppi_mmmm(vw$newY, rep(2000, nrow(vw$newY)), beta0 = m$beta0)
  est_direct <- ppi_mmmm(x, rep(2000, nrow(x)), beta0 = m$beta0, w = vw$w)
  expect_equal(est_direct, est_sim)
})

test_that("estimator1 and SE matches for simulated weights", {
  est_sim <- estimator1(vw$newY, acut = acut, incb = 1, beta = m$beta0, computeSE = TRUE)
  est_direct <- estimator1(m$sample, acut = acut, incb = 1, beta = m$beta0, w = vw$w, computeSE = TRUE)
  expect_equal(est_direct, est_sim)
})

test_that("estimatorall1 matches for simulated weights", {
  est_sim <- estimatorall1(vw$newY, acut = acut, w = rep(1, nrow(vw$newY)))
  est_direct <- estimatorall1(m$sample, acut = acut, w = vw$w)
  expect_equal(est_direct, est_sim)
})

test_that("dir_sqrt_prodh matches for simulated weights", {
  est_sim <- dir_sqrt_prodh(vw$newY, acut = acut, w = rep(1, nrow(vw$newY)))
  est_direct <- dir_sqrt_prodh(m$sample, acut = acut, w = vw$w)
  expect_equal(est_direct, est_sim)
})

test_that("ppi_alr_gengamma matches for simulated weights", {
  est_sim <- ppi_alr_gengamma(vw$newY, betap = m$beta0[3], w = rep(1, nrow(vw$newY)))
  est_direct <- ppi_alr_gengamma(m$sample, betap = m$beta0[3], w = vw$w)
  expect_equal(est_direct$est$paramvec, est_sim$est$paramvec)
})

test_that("dir_sqrt_minimah matches for simulated weights", {
  est_sim <- dir_sqrt_minimah(vw$newY, acut = acut)
  est_direct <- dir_sqrt_minimah(m$sample, acut = acut, w = vw$w)
  expect_equal(est_direct, est_sim)
})

test_that("dir_sqrt_prodh matches for simulated weights", {
  est_sim <- dir_sqrt_prodh(vw$newY, acut = acut)
  est_direct <- dir_sqrt_prodh(m$sample, acut = acut, w = vw$w)
  expect_equal(est_direct, est_sim)
})
