test_that("cppad_search gives similar result to cppad_closed", {
  m <- ppi_egmodel(100, maxden = 4)
  tapes <- buildsmotape("sphere", "ppi",
                        utape = rep(1/m$p, m$p),
                        usertheta = rep(NA, length(m$theta)),
                        weightname = "minsq", acut = 0.1)
  estsearch <- cppad_search(tapes$smotape, m$theta *0 + 1, m$sample, control = list(tol = 1E-12, maxit = 1000))
  estclosed <- cppad_closed(tapes$smotape, m$sample)
  expect_equal(estsearch$est, estclosed$est, tolerance = 1E-3)
  expect_equal(estsearch$SE, estclosed$SE, tolerance = 1E-3)
})

test_that("cppad_search with weights gives similar result to cppad_closed", {
  set.seed(124)
  m <- ppi_egmodel(100, maxden = 4)
  tapes <- buildsmotape("sphere", "ppi",
                        utape = rep(1/m$p, m$p),
                        usertheta = rep(NA, length(m$theta)),
                        weightname = "minsq", acut = 0.1)
  w <- runif(100)
  estsearch <- cppad_search(tapes$smotape, m$theta *0 + 1, m$sample, control = list(tol = 1E-12, maxit = 1000), w = w)
  estclosed <- cppad_closed(tapes$smotape, m$sample, w = w)
  expect_equal(estsearch$est, estclosed$est, tolerance = 1E-3)
  expect_equal(estsearch$SE, estclosed$SE)
})

test_that("cppad_search output value matches tape_smvalues result", {
  set.seed(1234)
  m <- ppi_egmodel(1000, maxden = 4)
  tapes <- buildsmotape("sphere", "ppi",
                        utape = rep(1/m$p, m$p),
                        usertheta = rep(NA, length(m$theta)),
                        weightname = "minsq", acut = 0.1)

  est <- cppad_search(tapes$smotape, m$theta *0 + 1, m$sample, control = list(tol = 1E-12, maxit = 1000))

  smvalues <- tape_smvalues_wsum(tapes$smotape, m$sample, est$est)
  expect_equal(est$value, smvalues$obj/nrow(m$sample))
  expect_equal(est$sqgradsize, sum(smvalues$grad^2))
})

