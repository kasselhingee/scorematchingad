m <- sec2_3model(10)

test_that("Correctly chooses Dirichlet", {
  out <- ppi(m$sample, AL = 0, bL = 0, man = "sphere", method = "direct",
             acut = 0.1,
             bdryweight = "minsq")
  expect_equal(out$fitfun, "estimator1_dir")

  out <- ppi(m$sample, AL = 0, bL = 0, man = "sphere", method = "direct",
             acut = 0.1,
             bdryweight = "prodsq")
  expect_equal(out$fitfun, "estimator2_dir")
})

test_that("Correctly chooses estimatorlog_ratio", {
  out <- ppi(m$sample, bL = 0, betap = tail(m$beta0, 1), man = "Ralr", method = "direct")
  expect_equal(out$fitfun, "estimatorlog_weight")
})

test_that("Correctly chooses sphere estimators with fixed beta", {
  out <- ppi(m$sample, bL = 0, beta = m$beta0, man = "sphere", method = "direct",
             acut = 0.1,
             bdryweight = "minsq")
  expect_equal(out$fitfun, "estimator1_zerob")
  out <- ppi(m$sample, bL = 0, beta = m$beta0, man = "sphere", method = "direct",
             acut = 0.1,
             bdryweight = "prodsq")
  expect_equal(out$fitfun, "estimator2_zerob")

  out <- ppi(m$sample, beta = m$beta0, man = "sphere", method = "direct",
             acut = 0.1,
             bdryweight = "minsq")
  expect_equal(out$fitfun, "estimator1_incb")
  out <- ppi(m$sample, beta = m$beta0, man = "sphere", method = "direct",
             acut = 0.1,
             bdryweight = "prodsq")
  expect_equal(out$fitfun, "estimator2_incb")
})

test_that("Correctly chooses sphere estimators for unfixed beta", {
  out <- ppi(m$sample, betap = tail(m$beta0, 1), man = "sphere", method = "direct",
             acut = 0.1,
             bdryweight = "minsq")
  expect_equal(out$fitfun, "estimatorall1_betap")
  out <- ppi(m$sample, man = "sphere", method = "direct",
             acut = 0.1,
             bdryweight = "minsq")
  expect_equal(out$fitfun, "estimatorall1_full")
})

test_that("Correctly chooses cppad", {
  # full direct with prodsq doesn't exist
  expect_warning(out <- ppi(m$sample, betap = tail(m$beta0, 1), man = "sphere", method = "direct",
             acut = 0.1,
             bdryweight = "prodsq"))
  expect_equal(out$fitfun, "cppad")

  expect_warning(out <- ppi(m$sample, man = "sphere", method = "direct",
             acut = 0.1,
             bdryweight = "prodsq"))

  out <- ppi(m$sample, beta = m$beta0, man = "sphere", method = "cppad",
             acut = 0.1,
             bdryweight = "prodsq")
  expect_equal(out$fitfun, "cppad")
})
