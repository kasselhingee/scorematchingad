# test of iterative solver for ppi (compare to cppad_closet()):
## microbiomfit with outliers

test_that("cppad_search goes towards cppad_closed result", {
  set.seed(354)
  m <- rppi_egmodel(100, maxden = 4)
  tapes <- tape_smd("sim","sqrt", "sph", "ppi",
                        ytape = rep(1/m$p, m$p),
                        usertheta = rep(NA, length(m$theta)),
                        bdryw = "minsq", acut = 0.1)
  estclosed <- cppad_closed(tapes$smdtape, m$sample)
  start <- (m$theta + estclosed$est)/2
  estsearch <- suppressWarnings(cppad_search(tapes$smdtape, start, m$sample, control = list(tol = 1E-5, maxit = 10)))
  expect_lt(sum((estclosed$est - estsearch$est)^2), sum((estclosed$est - start)^2))
})


