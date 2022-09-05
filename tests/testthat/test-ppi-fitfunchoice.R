m <- ppi_egmodel(10)

test_that("Correctly chooses Dirichlet", {
  out <- ppi(m$sample, ppi_paramvec(p=3, AL = 0, bL = 0), trans = "sqrt", method = "direct",
             acut = 0.1,
             bdryweight = "minsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "dir_sqrt_minimah")

  out <- ppi(m$sample, ppi_paramvec(p=3, AL = 0, bL = 0), trans = "sqrt", method = "direct",
             acut = 0.1,
             bdryweight = "prodsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "dir_sqrt_prodh")
})

test_that("Correctly chooses estimatorlog_ratio", {
  out <- ppi(m$sample, ppi_paramvec(p=3, bL = 0, betap = tail(m$beta0, 1)), trans = "alr", method = "direct")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "ppi_alr_gengamma")
})

test_that("Correctly chooses sphere estimators with fixed beta", {
  out <- ppi(m$sample, ppi_paramvec(bL = 0, beta = m$beta0), trans = "sqrt", method = "direct",
             acut = 0.1,
             bdryweight = "minsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "estimator1_zerob")
  out <- ppi(m$sample, ppi_paramvec(bL = 0, beta = m$beta0), trans = "sqrt", method = "direct",
             acut = 0.1,
             bdryweight = "prodsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "ppi_sqrt_prodh_zerob")

  out <- ppi(m$sample, ppi_paramvec(beta = m$beta0), trans = "sqrt", method = "direct",
             acut = 0.1,
             bdryweight = "minsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "estimator1_incb")
  out <- ppi(m$sample, ppi_paramvec(beta = m$beta0), trans = "sqrt", method = "direct",
             acut = 0.1,
             bdryweight = "prodsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "ppi_sqrt_prodh")
})

test_that("Correctly chooses sphere estimators for unfixed beta", {
  out <- ppi(m$sample, ppi_paramvec(p=3, betap = tail(m$beta0, 1)), trans = "sqrt", method = "direct",
             acut = 0.1,
             bdryweight = "minsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "ppi_sqrt_minimah_betap")
  out <- ppi(m$sample, trans = "sqrt", method = "direct",
             acut = 0.1,
             bdryweight = "minsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "ppi_sqrt_minimah_full")
})

test_that("Correctly chooses cppad", {
  # full direct with prodsq doesn't exist
  expect_warning(out <- ppi(m$sample, ppi_paramvec(p=3, betap = tail(m$beta0, 1)), trans = "sqrt", method = "direct",
             acut = 0.1,
             bdryweight = "prodsq"))
  expect_equal(out$info$method, "cppad")

  expect_warning(out <- ppi(m$sample, trans = "sqrt", method = "direct",
             acut = 0.1,
             bdryweight = "prodsq"))

  out <- ppi(m$sample, ppi_paramvec(beta = m$beta0), trans = "sqrt", method = "cppad",
             acut = 0.1,
             bdryweight = "prodsq")
  expect_equal(out$info$method, "cppad")
})
