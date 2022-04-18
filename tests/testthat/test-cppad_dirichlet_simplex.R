
test_that("cppad-based Score2 estimate leads to a match for large number of observations", {
  smofun <- ptapesmo_simplex(c(1,1,1,3,3,3), 3)
  beta = c(-0.3, -0.1, 3)
  n = 100
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  out <- optim(par = beta*0,
               fn = function(beta){smobj(smofun, beta,utabl)},
               gr = function(beta){smobjgrad(smofun, beta, utabl)},
               method = "BFGS")

  expect_equal(out$par, beta, tolerance = 1E-3, ignore_attr = TRUE)
})
