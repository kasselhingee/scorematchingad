test_that("estimatorlog_weight matches CppAD method for constant weight, p = 3", {
  set.seed(1234)
  m <- sec2_3model(1000, maxden = 4)

  est_cppad <- ppi_cppad(m$sample, bL = rep(0, 3-1), betap = m$beta0[3], man = "Ralr", weightname = "ones",
                         control = list(tol = 1E-10))

  est_direct <- estimatorlog_weight(m$sample, betap = m$beta0[3], weightW = rep(1, nrow(m$sample)))
  expect_equal(est_direct$ppi, est_cppad$est$theta, tolerance = 1E-5)
})

test_that("estimatorlog_weight matches CppAD method for constant weight, p = 5", {
  set.seed(1273)
  p = 5
  ALs <- exp(rsymmetricmatrix(p-1, -4, 4))
  bL <- rep(0, p-1)
  beta <- c(-0.7, -0.8, -0.3, 0, 0)
  prop <- rhybrid(1000, p, beta, ALs, bL, 35)$samp3

  est_cppad <- ppi_cppad(prop, bL = bL, betap = beta[p], man = "Ralr", weightname = "ones",
                         control = list(tol = 1E-20))
  expect_absdiff_lte_v(est_cppad$est$ALs, ALs, 3 * est_cppad$SE$ALs)

  est_direct <- estimatorlog_weight(prop, betap = beta[p], weightW = rep(1, nrow(prop)))
  # Get SE of this estimate using CppAD
  thetain <- ppi_cppad_thetaprocessor(p, bL = bL, betap = beta[p])
  tapes <- buildsmotape("Ralr", "ppi",
                        rep(1/p, p), thetain,
                        weightname = "ones",
                        acut = 1, verbose = FALSE)
  est_direct_SE <- smestSE(tapes$smotape, est_direct$ppi[is.na(thetain)], prop)
  expect_absdiff_lte_v(est_direct$ppi[is.na(thetain)], toPPIparamvec(ALs, bL, beta)[is.na(thetain)],
                       3 * est_direct_SE)

  smobjgrad(tapes$smotape, est_direct$ppi[is.na(thetain)], prop)
  smobjgrad(tapes$smotape, est_cppad$est$theta$ppi[is.na(thetain)], prop)


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
