test_that("ppi_alr_gengamma matches CppAD method for constant weight, p = 3", {
  set.seed(1234)
  m <- ppi_egmodel(100, maxden = 4)

  est_cppad <- ppi(m$sample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta0[3]), trans = "alr", method = "cppad", divweight = "ones",
                         control = list(tol = 1E-10))
  est_hardcoded <- ppi(m$sample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta0[3]), trans = "alr", method = "hardcoded")

  expect_equal(est_hardcoded$est$paramvec, est_cppad$est$paramvec)
})

test_that("Direct estimate has low smgrad values", {
  mnongamma <- ppi_egmodel(1)
  theta <- ppi_paramvec(beta = c(-0.95, -0.9, 0.5), AL = mnongamma$AL, bL = 0)
  set.seed(1234)
  Ycts <- rppi(1000, paramvec = theta)
  dsample <- round(Ycts * 100)/ 100
  dsample[, 3] <- 1 - rowSums(dsample[, 1:2])
  colMeans(dsample == 0)
  mean(apply(dsample, 1, min) == 0)  #0.96

  usertheta = ppi_paramvec(p = ncol(dsample),
                           bL = 0,
                           betap = tail(theta, 1))
  est_hardcoded <- ppi(dsample, paramvec = usertheta, trans = "alr", method = "hardcoded")

  hardcodedvals <- ppi_cppad_values(dsample,
         stheta = est_hardcoded$est$paramvec,
         isfixed = t_u2i(usertheta),
         man = "Ralr",
         hsqfun = "ones", 
         acut = 1)
  modelvals <- ppi_cppad_values(dsample,
         stheta = theta,
         isfixed = t_u2i(usertheta),
         man = "Ralr",
         hsqfun = "ones", 
         acut = 1)

  expect_lt(hardcodedvals$obj, modelvals$obj) #because for any given sample the estimate would be better than the true value

  expect_lt_v(abs(hardcodedvals$grad), abs(modelvals$grad)) #because the estimate will be better for any given sample
  expect_lt_v(abs(hardcodedvals$grad), rep(1E-10, length(hardcodedvals$grad)))
})

test_that("ppi_alr_gengamma matches CppAD method for constant weight and data with zeros, p = 3", {
  set.seed(1234)
  m <- ppi_egmodel(100, maxden = 4)
  dsample <- round(m$sample * 100)/100
  dsample[, 3] <- 1 - rowSums(dsample[, 1:2])

  est_hardcoded <- ppi(dsample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta0[3]), trans = "alr", method = "hardcoded")
  est_cppad <- ppi(dsample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta0[3]), trans = "alr", method = "cppad", divweight = "ones",
                         bdrythreshold = 1E-200,
                         control = list(tol = 1E-20))

  expect_equal(est_hardcoded$est$paramvec, est_cppad$est$paramvec)
})

test_that("ppi_alr_gengamma matches CppAD method for constant weight, p = 5", {
  skip_on_cran() #slow but useful
  set.seed(1273)
  p = 5
  ALs <- rsymmetricmatrix(p-1, -4, 4)
  bL <- rep(0, p-1)
  beta <- c(-0.7, -0.8, -0.3, 0, 0)
  # set.seed(1345) #this seed leads to samples with that give reasonable estimates
  set.seed(1111) #this seed leads to some ginormous elements for the second diagonal element of ALs
  prop <- rppi(1000, beta=beta, AL=ALs, bL=bL, maxden=5) #rppi_singly took 1005 seconds, rppi() took 13seconds

  est_cppad <- ppi(prop, ppi_paramvec(bL = bL, betap = beta[p]), trans = "alr", method = "cppad", divweight = "ones",
                         bdrythreshold = 1E-20,
                         control = list(tol = 1E-10))
  expect_absdiff_lte_v(est_cppad$est$ALs, ALs, 3 * est_cppad$SE$ALs)
  expect_absdiff_lte_v(est_cppad$est$beta, beta, 3 * est_cppad$SE$beta)
  #expect that the SE are small relative to size of the coefficients
  expect_lt(median(abs(est_cppad$SE$paramvec/est_cppad$est$paramvec), na.rm = TRUE), 0.3)

  est_hardcoded <- ppi_alr_gengamma(prop, betap = beta[p], w = rep(1, nrow(prop)))

  # Get SE of this estimate using CppAD
  thetain <- ppi_paramvec(p, bL = bL, betap = beta[p])
  tapes <- buildsmotape("Ralr", "ppi",
                        rep(1/p, p), thetain,
                        weightname = "ones",
                        acut = 1, verbose = FALSE)
  est_hardcoded_SE <- cppadSE(tapes$smotape, est_hardcoded$est$paramvec[is.na(thetain)], prop)
  expect_absdiff_lte_v(est_hardcoded$est$paramvec[is.na(thetain)], toPPIparamvec(ALs, bL, beta)[is.na(thetain)],
                       3 * est_hardcoded_SE)

  # check that hardcoded estimates are good according to smval and smvalgrad
  expect_lt(sum(smobjgrad(tapes$smotape, est_hardcoded$est$paramvec[is.na(thetain)], prop)^2), 1E-20)
  expect_lt(smobj(tapes$smotape, est_hardcoded$est$paramvec[is.na(thetain)], prop), est_cppad$info$smval)

  # check that estimates via cppad are close to hardcoded
  expect_absdiff_lte_v(est_hardcoded$est$paramvec[is.na(thetain)], est_cppad$est$paramvec[is.na(thetain)],
                       1.2 * est_hardcoded_SE) #the smovals are quite flat in the ALs dimensions for this region!
  # and that the beta estimates are really close to each other
  expect_equal(fromPPIparamvec(est_hardcoded$est$paramvec)$beta, est_cppad$est$beta, tolerance = 1E-3)
})


test_that("ppi_alr_gengamma matches for simulated weights", {
  set.seed(1234)
  m <- ppi_egmodel(100, maxden = 4)
  #simulate weights
  ind <- sample(1:100, 150, replace = TRUE)
  weights <- rep(0, 100)
  weights[as.numeric(names(table(ind)))] <- table(ind)
  newsample <- m$sample[ind, ]

  est_sim <- ppi_alr_gengamma(newsample, betap = m$beta0[3], w = rep(1, nrow(newsample)))
  est_hardcoded <- ppi_alr_gengamma(m$sample, betap = m$beta0[3], w = weights)
  expect_equal(est_hardcoded$est$paramvec, est_sim$est$paramvec)
})

test_that("ppi_alr_gengamma() and cppad match for a randomly selected weight vector", {
  set.seed(1234)
  m <- ppi_egmodel(100, maxden = 4)
  set.seed(1212)
  w <- runif(100)

  est_cppad <- ppi(m$sample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta0[3]), trans = "alr", method = "cppad", divweight = "ones",
                         w = w,
                         control = list(tol = 1E-10))
  est_hardcoded <- ppi(m$sample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta0[3]), trans = "alr", method = "hardcoded", w = w)

  expect_equal(est_hardcoded$est$paramvec, est_cppad$est$paramvec)
})


