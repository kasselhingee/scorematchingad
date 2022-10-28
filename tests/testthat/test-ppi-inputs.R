test_that("Inputs to ppi() are processed into the correct theta", {
  p = 3
  thetalength <- ppithetalength(p)
  expect_equal(ppi_paramvec(p), rep(NA, thetalength))
  # AL checks
  expect_equal(ppi_paramvec(p, AL = 3),
               c(rep(3, p-1), rep(3, (p-2) * (p-1)/2), rep(NA, p-1 + p)))
  expect_equal(ppi_paramvec(p, AL = NA), rep(NA, thetalength))
  expect_equal(ppi_paramvec(p, AL = "diag"),
               c(rep(NA, p-1), rep(0, (p-2) * (p-1)/2), rep(NA, p-1 + p)))
  expect_equal(ppi_paramvec(p, AL =  matrix(c(-166, 117, 117, -333), ncol = 2, nrow = 2)),
               c(-166, -333, 117, rep(NA, p-1 + p)))

  # A check
  Qin <- orthogmatwith111vec()
  Astar <- Qin %*% diag(c(-100, -200, 0)) %*% t(Qin)
  expect_equal(ppi_paramvec(p, Astar = Astar)[(thetalength - p + 1):thetalength], rep(NA_real_, p))
  expect_warning(expect_error(ppi_paramvec(p, AL = 0, Astar = 0)))

  # bL check
  expect_equal(ppi_paramvec(p, bL = NA),
               rep(NA, thetalength))
  expect_equal(ppi_paramvec(p, bL = 1),
               c(rep(NA, p-1 + (p-2) * (p-1)/2), rep(1, p-1), rep(NA,p)))
  expect_equal(ppi_paramvec(p, bL = c(2, 3)),
               c(rep(NA, p-1 + (p-2) * (p-1)/2), 2, 3, rep(NA,p)))
  expect_error(ppi_paramvec(p, bL = c(2, 3, 4)))

  # betaL check
  expect_equal(ppi_paramvec(p, betaL = NA),
               rep(NA, thetalength))
  expect_equal(ppi_paramvec(p, betaL = 1),
               c(rep(NA, p-1 + (p-2) * (p-1)/2), rep(NA, p-1), rep(1,p - 1), NA))
  expect_equal(ppi_paramvec(p, betaL = c(2, 3)),
               c(rep(NA, p-1 + (p-2) * (p-1)/2), rep(NA, p-1), 2,3, NA))
  expect_error(ppi_paramvec(p, betaL = c(2, 3, 4)))

  #betap check
  expect_equal(ppi_paramvec(p, betap = NA),
               rep(NA, thetalength))
  expect_equal(ppi_paramvec(p, betap = 1),
               c(rep(NA, thetalength - 1), 1))
  expect_error(ppi_paramvec(p, betap = c(1, 2)))

  #beta by itself check
  expect_equal(ppi_paramvec(p, beta = NA),
               rep(NA, thetalength))
  expect_equal(ppi_paramvec(p, beta = c(1.0, 2, NA)),
               c(rep(NA, thetalength - p), c(1.0,2,NA)))
  expect_error(ppi_paramvec(p, beta = c(1, 2)))
  expect_error(ppi_paramvec(p, beta = NA, betap = 1))
  expect_error(ppi_paramvec(p, beta = NA, betaL = c(1, 2)))

  #check supplying all together
  expect_equal(ppi_paramvec(p,
                           AL =  matrix(c(-166, 117, 117, -333), ncol = 2, nrow = 2),
                           bL = 2,
                           betaL = c(2, 3),
                           betap = 10),
    c(-166, -333, 117, rep(2, p-1), 2,3, 10))

  #p = 4
  p = 4
  expect_equal(ppi_paramvec(p, AL = "diag", bL = 0, betap = -0.5),
               c(rep(NA, p-1), rep(0, (p-2) * (p-1)/2), rep(0, p-1), rep(NA, p-1), -0.5))
})

test_that("ppi with cppad method works easily on ppi_egmodel", {
  set.seed(1245)
  model <- ppi_egmodel(100)
  out <- ppi(model$sample, trans = "sqrt", bdryweight = "minsq", acut = 0.1, method = "cppad", control = list(tol = 1E-10))
  expect_absdiff_lte_v(out$est$paramvec, model$theta, 3 * out$SE$paramvec)

  # try fixing betap
  out <- ppi(model$sample, ppi_paramvec(p = 3, betap = -0.5), trans = "sqrt", bdryweight = "minsq", acut = 0.1, method = "cppad", control = list(tol = 1E-10))
  expect_lte_v(abs(out$est$paramvec - model$theta), 3 * out$SE$paramvec)
  expect_equal(out$est$beta[model$p], -0.5)
  expect_equal(out$SE$beta[model$p], 0)

  # try fixing AL to diagonal
  AL = diag(c(-100, -50))
  beta = c(-0.8, -0.8, -0.5)
  bL = rep(0, 2)
  prop <- rppi(100, beta, AL, bL, maxden=4)
  theta <- toPPIparamvec(AL, bL, beta)
  out <- ppi(prop, ppi_paramvec(p=3, AL = "diag", betap = -0.5), trans = "sqrt", bdryweight = "minsq", acut = 0.1, method = "cppad", control = list(tol = 1E-10))
  expect_lte_v(abs(out$est$paramvec - theta), 3 * out$SE$paramvec)
  expect_equal(out$est$beta[model$p], -0.5)
  expect_equal(out$est$ALs[1, 2], 0)

  # try fixing AL to diagonal, bL to 0, betap = -0.5, on Ralr
  AL = diag(c(-100, -50))
  beta = c(-0.8, -0.8, -0.5)
  bL = rep(0, 2)
  prop <- rppi(100, beta, AL, bL, maxden=4)
  theta <- toPPIparamvec(AL, bL, beta)
  out <- ppi(prop, ppi_paramvec(p=3, AL = "diag", bL = 0, betap = -0.5), trans = "alr", bdryweight = "ones", method = "cppad", control = list(tol = 1E-10))
  expect_lte_v(abs(out$est$paramvec - theta), 3 * out$SE$paramvec)
  expect_equal(out$est$beta[model$p], -0.5)
  expect_equal(out$est$ALs[1, 2], 0)
})

test_that("ppi() uses paramvec_start}", {
  set.seed(1245)
  model <- ppi_egmodel(100)
  direct <- ppi(model$sample, trans = "sqrt", bdryweight = "minsq", acut = 0.1, method = "direct")
  out <- ppi(model$sample, trans = "sqrt", bdryweight = "minsq", acut = 0.1, method = "cppad", control = list(tol = 1E-10), paramvec_start = direct$est$paramvec)

  #expect very few iterations
  expect_lte_v(out$info$counts, rep(1, 1))
})

test_that("paramvec and paramvec_start are tested and made consistent correctly", {
  #basic elements to plug in
  AL <-matrix(c(-166, 117, 117, -333), ncol = 2, nrow = 2)
  beta <- c(0.5, -0.1, 0)
  bL <- rep(0, 2)
  p <- 3

  paramvec <- ppi_paramvec(AL = AL)
  paramvec_start <- ppi_paramvec(beta = beta)
  expect_error(t_us2s(paramvec, paramvec_start), regexp = "paramvec_start needs to supply.*4.*5")

  paramvec <- ppi_paramvec(AL = AL)
  paramvec_start <- ppi_paramvec(bL = bL, beta = beta)
  expect_equal(t_us2s(paramvec, paramvec_start), ppi_paramvec(AL = AL, bL = bL, beta = beta))

  paramvec <- ppi_paramvec(AL = AL)
  paramvec_start <- ppi_paramvec(AL = AL + 1, bL = bL, beta = beta)
  expect_warning(t_us2s(paramvec, paramvec_start), regexp = "paramvec_start inconsistent.*1.*2.*3")

  paramvec <- ppi_paramvec(AL = AL, bL = bL)
  paramvec_start <- ppi_paramvec(AL = AL, bL = bL + 1, beta = beta)
  expect_warning(t_us2s(paramvec, paramvec_start), regexp = "paramvec_start inconsistent.*4.*5")

  paramvec <- ppi_paramvec(AL = AL, beta = beta)
  paramvec_start <- ppi_paramvec(AL = AL, bL = bL, beta = beta + 1)
  expect_warning(t_us2s(paramvec, paramvec_start), regexp = "paramvec_start inconsistent.*6.*7.*8")
  expect_equal(suppressWarnings(t_us2s(paramvec, paramvec_start)), ppi_paramvec(AL = AL, bL = bL, beta = beta))
})
