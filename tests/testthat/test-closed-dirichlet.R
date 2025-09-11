skip_on_cran()

test_that("prodsq weights match estimator2", {
  acut = 0.1
  p = 3

  #simulate
  beta = c(-0.3, -0.1, 3)
  n = 10
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  tapes <- tape_smi(manifold = "sph",
                    uld = tape_uld_inbuilt("dirichlet", amdim = p),
                    transform = "sqrt",
                    bdryw = tape_bdryw_inbuilt("prodsq", rep(1/p, p), acut = acut))
  out <- cppad_closed(tapes$smi, Y = utabl)

  hardcodedestimate <- dir_sqrt_prodh(utabl, acut)
  expect_equal(out$est, hardcodedestimate, ignore_attr = TRUE)
})



test_that("minsq weights match estimator2", {
  acut = 0.1
  p = 3
  tapes <- tape_smi(manifold = "sph",
                    uld = tape_uld_inbuilt("dirichlet", amdim = p),
                    transform = "sqrt",
                    bdryw = tape_bdryw_inbuilt("minsq", rep(1/p, p), acut = acut))

  beta = c(-0.3, -0.1, 3)
  n = 10
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)
  out <- cppad_closed(tapes$smi, Y = utabl)

  hardcodedestimate <- dir_sqrt_minimah(utabl, acut)
  expect_equal(out$est, hardcodedestimate, ignore_attr = TRUE)
})


test_that("minsq weights match estimator2 for d = 4", {
  acut = 0.1
  p = 4
  tapes <- tape_smi(manifold = "sph",
                    uld = tape_uld_inbuilt("dirichlet", amdim = p),
                    transform = "sqrt",
                    bdryw = tape_bdryw_inbuilt("minsq", rep(1/p, p), acut = acut))

  beta = c(-0.3, -0.1, -0.2, 3)
  n = 10
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  out <- cppad_closed(tapes$smi, Y = utabl)
  hardcodedestimate <- dir_sqrt_minimah(utabl, acut)
  expect_equal(out$est, hardcodedestimate, ignore_attr = TRUE)
})

test_that("fixed beta[p] with minsq weights match true value", {
  acut = 0.1
  beta = c(-0.2, -0.1, 3)
  n = 10000
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  p = length(beta)
  tapes <- tape_smi(manifold = "sph",
                    uld = tape_uld_inbuilt("dirichlet", amdim = p),
                    transform = "sqrt",
                    fixedparams = c(NA, NA, beta[3]),
                    bdryw = tape_bdryw_inbuilt("minsq", rep(1/p, p), acut = acut))
  out <- cppad_closed(tapes$smi, Y = utabl)

  expect_absdiff_lte_v(out$est, beta[-p], 3 * out$SE)
})


test_that("cppad-based Score2 estimate leads to a match for large number of observations", {
  p = 3
  tapes <- tape_smi(manifold = "sim",
                    uld = tape_uld_inbuilt("dirichlet", amdim = p),
                    bdryw = tape_bdryw_inbuilt("prodsq", rep(1/p, p), acut = 0.1))

  beta = c(-0.3, -0.1, 3)
  n = 1000
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)
  out <- cppad_closed(tapes$smi, Y = utabl)
  expect_absdiff_lte_v(out$est, beta, out$SE * 3)
})

test_that("Simplex calculations are historically consistent", {
  p = 3
  tapes <- tape_smi(manifold = "sim",
                    uld = tape_uld_inbuilt("dirichlet", amdim = p),
                    bdryw = tape_bdryw_inbuilt("prodsq", rep(1/p, p), acut = 1))

  beta = c(-0.3, -0.1, 3)
  n = 10
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  smvals <- smvalues_wsum(tapes$smi, utabl, beta+0.5)
  expect_snapshot_value(smvals$obj/n, style = "json2", tolerance = 1E-5)
  expect_snapshot_value(smvals$grad/n, style = "json2", tolerance = 1E-5)
})
