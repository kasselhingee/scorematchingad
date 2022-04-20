
test_that("ADFun XPtr for computing values", {
  smofun <- ptapesmo(c(1,1,1,3,3,3), 3, manifoldname = "sphere", "prodsq") #tape of the score function
  beta = c(-0.3, -0.1, 3)
  n = 10
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- optim(par = beta*0,
               fn = function(beta){smobj(smofun, beta, utabl)},
               gr = function(beta){smobjgrad(smofun, beta, utabl)},
               method = "BFGS")

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimator2_dir(utabl, 0.001)
  expect_equal(out$par, directestimate, tolerance = 1E-3, ignore_attr = TRUE)
})

