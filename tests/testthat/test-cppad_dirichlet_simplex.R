
test_that("cppad-based Score2 estimate leads to a match for large number of observations", {
  smofun <- ptapesmo(c(1,1,1,3,3,3), 3,
                     manifoldname = "simplex", weightname = "prodsq",
                     acut = 0.01)
  beta = c(-0.3, -0.1, 3)
  n = 1000
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  out <- optim(par = beta*0,
               fn = function(beta){smobj(smofun, beta,utabl)},
               gr = function(beta){smobjgrad(smofun, beta, utabl)},
               method = "BFGS")

  expect_equal(out$par, beta, tolerance = 1E-1, ignore_attr = TRUE)
})

test_that("Simplex calculations are historically consistent", {
  smofun <- ptapesmo(c(1,1,1,3,3,3), 3,
                     manifoldname = "simplex", weightname = "prodsq",
                     acut = 1)
  beta = c(-0.3, -0.1, 3)
  n = 10
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  expect_snapshot_value(smobj(smofun, beta + 0.5, utabl), style = "json2", tolerance = 1E-5)
  expect_snapshot_value(smobjgrad(smofun, beta + 0.5, utabl), style = "json2", tolerance = 1E-5)
})
