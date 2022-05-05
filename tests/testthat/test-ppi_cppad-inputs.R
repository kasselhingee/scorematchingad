test_that("Inputs to ppi_cppad are processed into the correct theta", {
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
  expect_error(ppi_cppad_thetaprocessor(p, A = 0))
  expect_warning(expect_error(ppi_cppad_thetaprocessor(p, AL = 0, A = 0)))

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

  #check supplying all together
  expect_equal(ppi_cppad_thetaprocessor(p,
                           AL =  matrix(c(-166, 117, 117, -333), ncol = 2, nrow = 2),
                           bL = 2,
                           betaL = c(2, 3),
                           betap = 10),
    c(-166, -333, 117, rep(2, p-1), 2,3, 10))
})


