test_that("estimatorlog_weight matches CppAD method for constant weight", {
  set.seed(1234)
  m <- sec2_3model(1000, maxden = 4)

  est_cppad <- ppi_cppad(m$sample, bL = rep(0, 3-1), betap = m$beta0[3], man = "Ralr", weightname = "ones",
                         control = list(tol = 1E-10))

  est_direct <- estimatorlog_weight(m$sample, betap = m$beta0[3], weightW = rep(1, nrow(m$sample)))
  expect_equal(est_direct$ppi, est_cppad$est$theta, tolerance = 1E-5)
})

test_that("estimatorlog_weight matches for simulated weights", {
  set.seed(1234)
  m <- sec2_3model(1000, maxden = 4)
  #simulate weights
  ind <- sample(1:1000, 1500, replace = TRUE)
  weights <- rep(0, 1000)
  weights[as.numeric(names(table(ind)))] <- table(ind)
  newsample <- m$sample[ind, ]

  est_sim <- estimatorlog_weight(newsample, betap = m$beta0[3], weightW = rep(1, nrow(newsample)))
  est_direct <- estimatorlog_weight(m$sample, betap = m$beta0[3], weightW = weights)
  expect_equal(est_direct$ppi, est_sim$ppi)
})
