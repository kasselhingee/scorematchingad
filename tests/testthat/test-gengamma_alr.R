test_that("estimatorlog_weight matches CppAD method for constant weight, p = 3", {
  skip_on_cran() #v slow
  set.seed(1234)
  m <- sec2_3model(10000, maxden = 4)

  est_cppad <- ppi(m$sample, bL = rep(0, 3-1), betap = m$beta0[3], man = "Ralr", method = "cppad", bdryweight = "ones",
                         control = list(tol = 1E-10))

  est_direct <- estimatorlog_weight(m$sample, betap = m$beta0[3], weightW = rep(1, nrow(m$sample)))
  expect_equal(est_direct$ppi, est_cppad$est$theta, tolerance = 1E-5)
})

test_that("estimatorlog_weight matches CppAD method for constant weight, p = 5", {
  skip_on_cran() #slow but useful
  set.seed(1273)
  p = 5
  ALs <- rsymmetricmatrix(p-1, -4, 4)
  bL <- rep(0, p-1)
  beta <- c(-0.7, -0.8, -0.3, 0, 0)
  # set.seed(1345) #this seed leads to samples with that give reasonable estimates
  set.seed(1111) #this seed leads to some ginormous elements for the second diagonal element of ALs
  prop <- rppi(1000, p, beta, ALs, bL, 5)$samp3 #rppi_singly took 1005 seconds, rppi() took 13seconds

  est_cppad <- ppi(prop, bL = bL, betap = beta[p], man = "Ralr", method = "cppad", bdryweight = "ones",
                         bdrythreshold = 1E-20,
                         control = list(tol = 1E-10))
  expect_absdiff_lte_v(est_cppad$est$ALs, ALs, 3 * est_cppad$SE$ALs)
  expect_absdiff_lte_v(est_cppad$est$beta, beta, 3 * est_cppad$SE$beta)
  #expect that the SE are small relative to size of the coefficients
  expect_lt(median(abs(est_cppad$SE$theta/est_cppad$est$theta), na.rm = TRUE), 0.3)

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

  # check that direct estimates are good according to smval and smvalgrad
  expect_lt(sum(smobjgrad(tapes$smotape, est_direct$ppi[is.na(thetain)], prop)^2), 1E-20)
  expect_lt(smobj(tapes$smotape, est_direct$ppi[is.na(thetain)], prop), est_cppad$smval)

  # check that estimates via cppad are close to direct
  expect_absdiff_lte_v(est_direct$ppi[is.na(thetain)], est_cppad$est$theta[is.na(thetain)],
                       1.2 * est_direct_SE) #the smovals are quite flat in the ALs dimensions for this region!
  # and that the beta estimates are really close to each other
  expect_equal(fromPPIparamvec(est_direct$ppi)$beta, est_cppad$est$beta, tolerance = 1E-3)
})


test_that("estimatorlog_weight matches for simulated weights", {
  set.seed(1234)
  m <- sec2_3model(100, maxden = 4)
  #simulate weights
  ind <- sample(1:100, 150, replace = TRUE)
  weights <- rep(0, 100)
  weights[as.numeric(names(table(ind)))] <- table(ind)
  newsample <- m$sample[ind, ]

  est_sim <- estimatorlog_weight(newsample, betap = m$beta0[3], weightW = rep(1, nrow(newsample)))
  est_direct <- estimatorlog_weight(m$sample, betap = m$beta0[3], weightW = weights)
  expect_equal(est_direct$ppi, est_sim$ppi)
})
