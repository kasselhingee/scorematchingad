test_that("Inputs to ppi() are processed into the correct theta", {
  p = 3
  thetalength <- ppithetalength(p)
  expect_equal(ppi_cppad_thetaprocessor(p), rep(NA, thetalength))
  # AL checks
  expect_equal(ppi_cppad_thetaprocessor(p, AL = 3),
               c(rep(3, p-1), rep(3, (p-2) * (p-1)/2), rep(NA, p-1 + p)))
  expect_equal(ppi_cppad_thetaprocessor(p, AL = NA), rep(NA, thetalength))
  expect_equal(ppi_cppad_thetaprocessor(p, AL = "diag"),
               c(rep(NA, p-1), rep(0, (p-2) * (p-1)/2), rep(NA, p-1 + p)))
  expect_equal(ppi_cppad_thetaprocessor(p, AL =  matrix(c(-166, 117, 117, -333), ncol = 2, nrow = 2)),
               c(-166, -333, 117, rep(NA, p-1 + p)))

  # A check
  Qin <- orthogmatwith111vec()
  Astar <- Qin %*% diag(c(-100, -200, 0)) %*% t(Qin)
  expect_equal(ppi_cppad_thetaprocessor(p, Astar = Astar)[(thetalength - p + 1):thetalength], rep(NA_real_, p))
  expect_warning(expect_error(ppi_cppad_thetaprocessor(p, AL = 0, Astar = 0)))

  # bL check
  expect_equal(ppi_cppad_thetaprocessor(p, bL = NA),
               rep(NA, thetalength))
  expect_equal(ppi_cppad_thetaprocessor(p, bL = 1),
               c(rep(NA, p-1 + (p-2) * (p-1)/2), rep(1, p-1), rep(NA,p)))
  expect_equal(ppi_cppad_thetaprocessor(p, bL = c(2, 3)),
               c(rep(NA, p-1 + (p-2) * (p-1)/2), 2, 3, rep(NA,p)))
  expect_error(ppi_cppad_thetaprocessor(p, bL = c(2, 3, 4)))

  # betaL check
  expect_equal(ppi_cppad_thetaprocessor(p, betaL = NA),
               rep(NA, thetalength))
  expect_equal(ppi_cppad_thetaprocessor(p, betaL = 1),
               c(rep(NA, p-1 + (p-2) * (p-1)/2), rep(NA, p-1), rep(1,p - 1), NA))
  expect_equal(ppi_cppad_thetaprocessor(p, betaL = c(2, 3)),
               c(rep(NA, p-1 + (p-2) * (p-1)/2), rep(NA, p-1), 2,3, NA))
  expect_error(ppi_cppad_thetaprocessor(p, betaL = c(2, 3, 4)))

  #betap check
  expect_equal(ppi_cppad_thetaprocessor(p, betap = NA),
               rep(NA, thetalength))
  expect_equal(ppi_cppad_thetaprocessor(p, betap = 1),
               c(rep(NA, thetalength - 1), 1))
  expect_error(ppi_cppad_thetaprocessor(p, betap = c(1, 2)))

  #beta by itself check
  expect_equal(ppi_cppad_thetaprocessor(p, beta = NA),
               rep(NA, thetalength))
  expect_equal(ppi_cppad_thetaprocessor(p, beta = c(1.0, 2, NA)),
               c(rep(NA, thetalength - p), c(1.0,2,NA)))
  expect_error(ppi_cppad_thetaprocessor(p, beta = c(1, 2)))
  expect_error(ppi_cppad_thetaprocessor(p, beta = NA, betap = 1))
  expect_error(ppi_cppad_thetaprocessor(p, beta = NA, betaL = c(1, 2)))

  #check supplying all together
  expect_equal(ppi_cppad_thetaprocessor(p,
                           AL =  matrix(c(-166, 117, 117, -333), ncol = 2, nrow = 2),
                           bL = 2,
                           betaL = c(2, 3),
                           betap = 10),
    c(-166, -333, 117, rep(2, p-1), 2,3, 10))

  #p = 4
  p = 4
  expect_equal(ppi_cppad_thetaprocessor(p, AL = "diag", bL = 0, betap = -0.5),
               c(rep(NA, p-1), rep(0, (p-2) * (p-1)/2), rep(0, p-1), rep(NA, p-1), -0.5))
})

test_that("ppi with cppad method works easily on sec2_3model", {
  set.seed(1245)
  model <- sec2_3model(1000)
  out <- ppi(model$sample, man = "sphere", bdryweight = "minsq", acut = 0.1, method = "cppad")
  cdabyppi::expect_lt_v(abs(out$est$theta - model$theta) / out$SE$theta, 3)

  # try fixing betap
  out <- ppi(model$sample, betap = -0.5, man = "sphere", bdryweight = "minsq", acut = 0.1, method = "cppad")
  cdabyppi:::expect_lte_v(abs(out$est$theta - model$theta), 3 * out$SE$theta)
  expect_equal(out$est$beta[model$p], -0.5)
  expect_equal(out$SE$beta[model$p], 0)

  # try fixing AL to diagonal
  AL = diag(c(-100, -50))
  beta = c(-0.8, -0.8, -0.5)
  bL = rep(0, 2)
  prop <- rhybrid(100, 3, beta, AL, bL, 4)[[1]]
  theta <- toPPIparamvec(AL, bL, beta)
  out <- ppi(prop, AL = "diag", betap = -0.5, man = "sphere", bdryweight = "minsq", acut = 0.1, method = "cppad")
  cdabyppi:::expect_lte_v(abs(out$est$theta - theta), 3 * out$SE$theta)
  expect_equal(out$est$beta[model$p], -0.5)
  expect_equal(out$est$ALs[1, 2], 0)

  # try fixing AL to diagonal, bL to 0, betap = -0.5, on Ralr
  AL = diag(c(-100, -50))
  beta = c(-0.8, -0.8, -0.5)
  bL = rep(0, 2)
  prop <- rhybrid(100, 3, beta, AL, bL, 4)[[1]]
  theta <- toPPIparamvec(AL, bL, beta)
  out <- ppi(prop, AL = "diag", bL = 0, betap = -0.5, man = "Ralr", bdryweight = "ones", method = "cppad")
  cdabyppi:::expect_lte_v(abs(out$est$theta - theta), 3 * out$SE$theta)
  expect_equal(out$est$beta[model$p], -0.5)
  expect_equal(out$est$ALs[1, 2], 0)
})
