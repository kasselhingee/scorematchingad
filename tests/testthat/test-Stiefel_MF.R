test_that("estimating matrix von Mises-Fisher get close to correct", {
  set.seed(1)
  M <- matrix(rnorm(4*2), nrow = 4, ncol = 2)
  set.seed(3)
  samp <- replicate(100, rstiefel::rmf.matrix(M))
  
  set.seed(4)
  Mother <- matrix(rnorm(nrow(M)*ncol(M)),
                   nrow = nrow(M), ncol = ncol(M))
  tapes <- tape_smi(
    manifold = make_manifold("Stiefel", param1 = nrow(M), param2 = ncol(M)),
    uld = "Stiefel_MF",
    xtape = vec(samp[,,1]),
    fixedparams = NA_real_*vec(M))
  expect_true(testquadratic(tapes$smi))
  Y <- t(apply(samp, 3, vec))
  est <- cppad_closed(tapes$smi, Y)
  expect_absdiff_lte_v(est$est, vec(M), 3 * est$SE)
})
