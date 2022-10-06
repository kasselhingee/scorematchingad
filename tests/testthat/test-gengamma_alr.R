test_that("ppi_alr_gengamma matches CppAD method for constant weight, p = 3", {
  set.seed(1234)
  m <- ppi_egmodel(100, maxden = 4)

  est_cppad <- ppi(m$sample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta0[3]), trans = "alr", method = "cppad", bdryweight = "ones",
                         control = list(tol = 1E-10))
  est_direct <- ppi(m$sample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta0[3]), trans = "alr", method = "direct")

  expect_equal(est_direct$est$paramvec, est_cppad$est$paramvec)
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
  est_direct <- ppi(dsample, paramvec = usertheta, trans = "alr", method = "direct")

  smotapes <- buildsmotape(manifoldname = "Ralr",
                         llname = "ppi",
                         utape = rep(1, ncol(dsample))/ncol(dsample),
                         usertheta = usertheta,
                         weightname = "ones", 
                         acut = 1,
                         thetatape_creator = function(n){seq(length.out = n)},
                         verbose = FALSE)
   
  eginteriorpt <- rep(1, ncol(dsample))/ncol(dsample) 
  theta_for_taping <- t_si2f(smotapes$info$starttheta, smotapes$info$isfixed)
  smofun_u <- swapDynamic(smotapes$smotape, eginteriorpt, theta_for_taping)
  Jsmofun_u <- pTapeJacobianSwap(smotapes$smotape, theta_for_taping, eginteriorpt)
  Hsmofun_u <- pTapeHessianSwap(smotapes$smotape, theta_for_taping, eginteriorpt)

  datasplit <- simplex_boundarysplit(dsample, 
                         bdrythreshold = 1E-10, 
                         shiftsize = 1E-10)

  objval <- smobj(smofun = smotapes$smotape,
        theta = t_si2f(est_direct$est$paramvec, smotapes$info$isfixed),
        utabl = datasplit$interior,
        smofun_u = smofun_u,
        uboundary = datasplit$uboundary,
        boundaryapprox = datasplit$boundaryapprox,
        approxorder = 10
        )

  objval_model <- smobj(smofun = smotapes$smotape,
        theta = t_si2f(theta, smotapes$info$isfixed),
        utabl = datasplit$interior,
        smofun_u = smofun_u,
        uboundary = datasplit$uboundary,
        boundaryapprox = datasplit$boundaryapprox,
        approxorder = 10
        )
  expect_lt(objval, objval_model) #because for any given sample the estimate would be better than the true value

  gradval <- smobjgrad(smofun = smotapes$smotape,
        theta = t_si2f(est_direct$est$paramvec, smotapes$info$isfixed),
        utabl = datasplit$interior,
        Jsmofun_u = Jsmofun_u,
        uboundary = datasplit$uboundary,
        boundaryapprox = datasplit$boundaryapprox,
        approxorder = 10
        )
  
  gradval_model <- smobjgrad(smofun = smotapes$smotape,
        theta = t_si2f(theta, smotapes$info$isfixed),
        utabl = datasplit$interior,
        Jsmofun_u = Jsmofun_u,
        uboundary = datasplit$uboundary,
        boundaryapprox = datasplit$boundaryapprox,
        approxorder = 10
        )
  expect_lt_v(abs(gradval), abs(gradval_model)) #because the estimate will be better for any given sample
  expect_lt_v(abs(gradval), rep(1E-10, length(gradval)))
})

test_that("model theta has low cppad smgrad values for rounded data with zeroes", {
  # my guess is that this works when the direct estimate is close to the model theta
  mnongamma <- ppi_egmodel(1)
  theta <- ppi_paramvec(beta = c(-0.95, -0.9, 0.5), AL = mnongamma$AL, bL = 0)
  set.seed(1234)
  Ycts <- rppi(10000, paramvec = theta)
  dsample <- round(Ycts * 100)/ 100
  dsample[, 3] <- 1 - rowSums(dsample[, 1:2])
  colMeans(dsample == 0)
  mean(apply(dsample, 1, min) == 0)

  usertheta = ppi_paramvec(p = ncol(dsample),
                           bL = 0,
                           betap = tail(theta, 1))
  est_direct <- ppi(dsample, paramvec = usertheta, trans = "alr", method = "direct")
  stopifnot(mean(abs(est_direct$est$paramvec - theta)) < 15)

  smotapes <- buildsmotape(manifoldname = "Ralr",
                         llname = "ppi",
                         utape = rep(1, ncol(dsample))/ncol(dsample),
                         usertheta = usertheta,
                         weightname = "ones", 
                         acut = 1,
                         thetatape_creator = function(n){seq(length.out = n)},
                         verbose = FALSE)
   
  eginteriorpt <- rep(1, ncol(dsample))/ncol(dsample) 
  theta_for_taping <- t_si2f(smotapes$info$starttheta, smotapes$info$isfixed)
  smofun_u <- swapDynamic(smotapes$smotape, eginteriorpt, theta_for_taping)
  Jsmofun_u <- pTapeJacobianSwap(smotapes$smotape, theta_for_taping, eginteriorpt)
  Hsmofun_u <- pTapeHessianSwap(smotapes$smotape, theta_for_taping, eginteriorpt)

  datasplit <- simplex_boundarysplit(dsample, 
                         bdrythreshold = 1E-10, 
                         shiftsize = 1E-10)

  objval <- smobj(smofun = smotapes$smotape,
        theta = t_si2f(theta, smotapes$info$isfixed),
        utabl = datasplit$interior,
        smofun_u = smofun_u,
        uboundary = datasplit$uboundary,
        boundaryapprox = datasplit$boundaryapprox,
        approxorder = 10
        )

  gradval <- smobjgrad(smofun = smotapes$smotape,
        theta = t_si2f(theta, smotapes$info$isfixed),
        utabl = datasplit$interior,
        Jsmofun_u = Jsmofun_u,
        uboundary = datasplit$uboundary,
        boundaryapprox = datasplit$boundaryapprox,
        approxorder = 10
        )

  expect_lt_v(abs(gradval), rep(1E-5, sum(!smotapes$info$isfixed)))

  gradval_direct_est <- smobjgrad(smofun = smotapes$smotape,
        theta = t_si2f(est_direct$est$paramvec, smotapes$info$isfixed),
        utabl = datasplit$interior,
        Jsmofun_u = Jsmofun_u,
        uboundary = datasplit$uboundary,
        boundaryapprox = datasplit$boundaryapprox,
        approxorder = 10
        )
  expect_lt_v(abs(gradval_direct_est), rep(1E-5, sum(!smotapes$info$isfixed)))
})

test_that("ppi_alr_gengamma matches CppAD method for constant weight and data with zeros, p = 3", {
  set.seed(1234)
  m <- ppi_egmodel(100, maxden = 4)
  dsample <- round(m$sample * 100)/100
  dsample[, 3] <- 1 - rowSums(dsample[, 1:2])

  est_direct <- ppi(dsample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta0[3]), trans = "alr", method = "direct")
  est_cppad <- ppi(dsample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta0[3]), trans = "alr", method = "cppad", bdryweight = "ones",
                         bdrythreshold = 1E-200,
                         control = list(tol = 1E-20))

  expect_equal(est_direct$est$paramvec, est_cppad$est$paramvec)
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

  est_cppad <- ppi(prop, ppi_paramvec(bL = bL, betap = beta[p]), trans = "alr", method = "cppad", bdryweight = "ones",
                         bdrythreshold = 1E-20,
                         control = list(tol = 1E-10), w = NULL) #w = NULL here to temporarily dodge the issue with weights of 1 generating different results to no weights
  expect_absdiff_lte_v(est_cppad$est$ALs, ALs, 3 * est_cppad$SE$ALs)
  expect_absdiff_lte_v(est_cppad$est$beta, beta, 3 * est_cppad$SE$beta)
  #expect that the SE are small relative to size of the coefficients
  expect_lt(median(abs(est_cppad$SE$paramvec/est_cppad$est$paramvec), na.rm = TRUE), 0.3)

  est_direct <- ppi_alr_gengamma(prop, betap = beta[p], w = rep(1, nrow(prop)))

  # Get SE of this estimate using CppAD
  thetain <- ppi_paramvec(p, bL = bL, betap = beta[p])
  tapes <- buildsmotape("Ralr", "ppi",
                        rep(1/p, p), thetain,
                        weightname = "ones",
                        acut = 1, verbose = FALSE)
  est_direct_SE <- cppadSE(tapes$smotape, est_direct$est$paramvec[is.na(thetain)], prop)
  expect_absdiff_lte_v(est_direct$est$paramvec[is.na(thetain)], toPPIparamvec(ALs, bL, beta)[is.na(thetain)],
                       3 * est_direct_SE)

  # check that direct estimates are good according to smval and smvalgrad
  expect_lt(sum(smobjgrad(tapes$smotape, est_direct$est$paramvec[is.na(thetain)], prop)^2), 1E-20)
  expect_lt(smobj(tapes$smotape, est_direct$est$paramvec[is.na(thetain)], prop), est_cppad$info$smval)

  # check that estimates via cppad are close to direct
  expect_absdiff_lte_v(est_direct$est$paramvec[is.na(thetain)], est_cppad$est$paramvec[is.na(thetain)],
                       1.2 * est_direct_SE) #the smovals are quite flat in the ALs dimensions for this region!
  # and that the beta estimates are really close to each other
  expect_equal(fromPPIparamvec(est_direct$est$paramvec)$beta, est_cppad$est$beta, tolerance = 1E-3)
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
  est_direct <- ppi_alr_gengamma(m$sample, betap = m$beta0[3], w = weights)
  expect_equal(est_direct$est$paramvec, est_sim$est$paramvec)
})

test_that("ppi_alr_gengamma() and cppad match for a randomly selected weight vector", {
  set.seed(1234)
  m <- ppi_egmodel(100, maxden = 4)
  set.seed(1212)
  w <- runif(100)

  est_cppad <- ppi(m$sample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta0[3]), trans = "alr", method = "cppad", bdryweight = "ones",
                         w = w,
                         control = list(tol = 1E-10))
  est_direct <- ppi(m$sample, ppi_paramvec(bL = rep(0, 3-1), betap = m$beta0[3]), trans = "alr", method = "direct", w = w)

  expect_equal(est_direct$est$paramvec, est_cppad$est$paramvec)
})


