skip_on_cran()

test_that("prodsq weights match estimator2", {
  acut = 0.1
  p = 3
  tapes <- buildsmotape("sphere", "dirichlet",
               rep(1, p), rep(NA, p),
               "prodsq", acut = acut)

  #simulate
  beta = c(-0.3, -0.1, 3)
  n = 10
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2017/11/why-optim-is-out-of-date/)

  out <- Rcgmin::Rcgmin(par = beta*0,
               fn = function(beta){smobj(tapes$smotape, beta, utabl)},
               gr = function(beta){smobjgrad(tapes$smotape, beta, utabl)})

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- dir_sqrt_prodh(utabl, acut)
  expect_equal(out$par, directestimate, tolerance = 1E-3, ignore_attr = TRUE)
})



test_that("minsq weights match estimator2", {
  acut = 0.1
  p = 3
  tapes <- buildsmotape("sphere", "dirichlet",
                        rep(1, p), rep(NA, p),
                        "minsq", acut = acut)

  beta = c(-0.3, -0.1, 3)
  n = 10
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- Rcgmin::Rcgmin(par = beta*0,
               fn = function(beta){smobj(tapes$smotape, beta, utabl)},
               gr = function(beta){smobjgrad(tapes$smotape, beta, utabl)})

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- dir_sqrt_minimah(utabl, acut)
  expect_equal(out$par, directestimate, tolerance = 1E-3, ignore_attr = TRUE)
})


test_that("minsq weights match estimator2 for d = 4", {
  acut = 0.1
  p = 4
  tapes <- buildsmotape("sphere", "dirichlet",
                        rep(1, p), rep(NA, p),
                        "minsq", acut = acut)

  beta = c(-0.3, -0.1, -0.2, 3)
  n = 10
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- Rcgmin::Rcgmin(par = beta*0,
               fn = function(beta){smobj(tapes$smotape, beta, utabl)},
               gr = function(beta){smobjgrad(tapes$smotape, beta, utabl)})

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- dir_sqrt_minimah(utabl, acut)
  expect_equal(out$par, directestimate, tolerance = 1E-3, ignore_attr = TRUE)
})

test_that("fixed beta[p] with minsq weights match true value", {
  acut = 0.1
  beta = c(-0.2, -0.1, 3)
  n = 10000
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  p = length(beta)
  tapes <- buildsmotape("sphere", "dirichlet",
                        rep(0.1, p), c(NA, NA, beta[3]),
                        "minsq", acut = acut)

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2017/11/why-optim-is-out-of-date/)
  out <- Rcgmin::Rcgmin(par = beta[-3] * 0,
               fn = function(beta){smobj(tapes$smotape, beta, utabl)},
               gr = function(beta){smobjgrad(tapes$smotape, beta, utabl)})

  directestimate <- dir_sqrt_minimah(utabl, acut)
  expect_equal(out$par, directestimate[-3], tolerance = 0.05, ignore_attr = TRUE) #tolerance is usually relative
  expect_equal(out$par, beta[-3], tolerance = 0.1, ignore_attr = TRUE) #tolerance is usually relative
})


test_that("cppad-based Score2 estimate leads to a match for large number of observations", {
  p = 3
  tapes <- buildsmotape("simplex", "dirichlet",
                        rep(0.1, p), rep(NA, p),
                        "prodsq", acut = 0.1)

  beta = c(-0.3, -0.1, 3)
  n = 1000
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  out <- Rcgmin::Rcgmin(par = beta*0,
               fn = function(beta){smobj(tapes$smotape, beta,utabl)},
               gr = function(beta){smobjgrad(tapes$smotape, beta, utabl)})

  expect_equal(out$par, beta, tolerance = 1E-1, ignore_attr = TRUE)
})

test_that("Simplex calculations are historically consistent", {
  p = 3
  tapes <- buildsmotape("simplex", "dirichlet",
                        rep(0.1, p), rep(NA, p),
                        "prodsq", acut = 1)

  beta = c(-0.3, -0.1, 3)
  n = 10
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  expect_snapshot_value(smobj(tapes$smotape, beta + 0.5, utabl), style = "json2", tolerance = 1E-5)
  expect_snapshot_value(smobjgrad(tapes$smotape, beta + 0.5, utabl), style = "json2", tolerance = 1E-5)
})
