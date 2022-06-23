#helper function: given a sample, simulate integer weights, then return weights and a version where each measurement is replicated by the weight size
virtualweights <- function(Y, sizefactor = 1.5){
  ind <- sample(1:nrow(Y), ceiling(sizefactor*nrow(Y)), replace = TRUE)
  weights <- rep(0, nrow(Y))
  weights[as.numeric(names(table(ind)))] <- table(ind)
  newsample <- Y[ind, ]
  return(list(
    newY = newsample,
    w = weights
  ))
}


test_that("estimatorlog_weight matches for simulated weights", {
  set.seed(1234)
  m <- sec2_3model(1000, maxden = 4)
  #simulate weights
  set.seed(134)
  vw <- virtualweights(m$sample)

  est_sim <- estimatorlog_weight(vw$newY, betap = m$beta0[3], weightW = rep(1, nrow(vw$newY)))
  est_direct <- estimatorlog_weight(m$sample, betap = m$beta0[3], weightW = vw$w)
  expect_equal(est_direct$ppi, est_sim$ppi)
})

test_that("dirichmom matches for simulated weights", {
  set.seed(1234)
  m <- sec2_3model(1000, maxden = 4)
  #simulate weights
  set.seed(134)
  vw <- virtualweights(m$sample)

  est_sim <- dirichmom(vw$newY)
  est_direct <- dirichmom(m$sample, w = vw$w)
  expect_equal(est_direct, est_sim)
})

test_that("estimator1_dir matches for simulated weights", {
  set.seed(1234)
  m <- sec2_3model(1000, maxden = 4)
  #simulate weights
  set.seed(134)
  vw <- virtualweights(m$sample)

  acut = 0.1
  est_sim <- estimator1_dir(vw$newY, acut = acut)
  est_direct <- estimator1_dir(m$sample, acut = acut, w = vw$w)
  expect_equal(est_direct, est_sim)
})
