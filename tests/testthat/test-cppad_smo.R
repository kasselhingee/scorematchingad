
test_that("prodsq weights match estimator2", {
  acut = 0.01
  psphere <- pmanifold("sphere")
  pdir <- ptapell(c(1,1,1), c(3,3,3), llname = "dirichlet", psphere, fixedtheta = c(FALSE, FALSE, FALSE))
  smofun <- ptapesmo(c(1,1,1),
                     c(3,3,3),
                     pll = pdir,
                     pman = psphere, "prodsq", acut = acut) #tape of the score function
  beta = c(-0.3, -0.1, 3)
  n = 10
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2017/11/why-optim-is-out-of-date/)
  out <- optim(par = beta*0,
               fn = function(beta){smobj(smofun, beta, utabl)},
               gr = function(beta){smobjgrad(smofun, beta, utabl)},
               method = "BFGS")

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimator2_dir(utabl, acut)
  expect_equal(out$par, directestimate, tolerance = 1E-3, ignore_attr = TRUE)
})



test_that("minsq weights match estimator2", {
  acut = 0.1
  psphere <- pmanifold("sphere")
  pdir <- ptapell(c(1,1,1), c(3,3,3), llname = "dirichlet", psphere, fixedtheta = rep(FALSE, 3))
  smofun <- ptapesmo(c(1,1,1),
                     c(3,3,3),
                     pll = pdir,
                     pman = psphere, "minsq", acut = acut) #tape of the score function
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
  directestimate <- estimator1_dir(utabl, acut)
  expect_equal(out$par, directestimate, tolerance = 1E-3, ignore_attr = TRUE)
})


test_that("minsq weights match estimator2 for d = 4", {
  acut = 0.1
  psphere <- pmanifold("sphere")
  pdir <- ptapell(c(1,1,1,1), c(3,3,3,3), llname = "dirichlet", psphere, fixedtheta = rep(FALSE, 4))
  smofun <- ptapesmo(c(1,1,1,1),
                     c(3,3,3,3),
                     pll = pdir,
                     pman = psphere, "minsq", acut = acut) #tape of the score function
  beta = c(-0.3, -0.1, -0.2, 3)
  n = 10
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- optim(par = beta*0,
               fn = function(beta){smobj(smofun, beta, utabl)},
               gr = function(beta){smobjgrad(smofun, beta, utabl)},
               method = "BFGS")

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimator1_dir(utabl, acut)
  expect_equal(out$par, directestimate, tolerance = 1E-3, ignore_attr = TRUE)
})

test_that("fixed beta[p] with minsq weights match true value", {
  acut = 0.1
  beta = c(-0.2, -0.1, 3)
  n = 10000
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  psphere <- pmanifold("sphere")
  pdir <- ptapell(c(0.1,0.1,0.1), c(1,2,3), llname = "dirichlet", psphere, fixedtheta = c(FALSE, FALSE, TRUE))
  smofun <- ptapesmo(c(0.1,0.1,0.1),
                     c(1, 2),
                     pll = pdir,
                     pman = psphere, "minsq", acut = acut) #tape of the score function

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2017/11/why-optim-is-out-of-date/)
  out <- optim(par = beta[-3] * 0,
               fn = function(beta){smobj(smofun, beta, utabl)},
               gr = function(beta){smobjgrad(smofun, beta, utabl)},
               method = "BFGS")

  directestimate <- estimator1_dir(utabl, acut)
  expect_equal(out$par, directestimate[-3], tolerance = 0.05, ignore_attr = TRUE) #tolerance is usually relative
  expect_equal(out$par, beta[-3], tolerance = 0.1, ignore_attr = TRUE) #tolerance is usually relative
})
