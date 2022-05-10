test_that("Fitting ppi via alr transform with fixed beta gets close to true values", {
  set.seed(1234)
  model <- sec2_3model(1000, maxden = 4)

  acut = 0.1 #not needed for Ralr per se, but code still expects it

  pman <- pmanifold("Ralr")
  pppi <- ptapell(rep(0.1, model$p - 1), model$theta, llname = "ppi", pman,
                  fixedtheta = c(rep(FALSE, length(model$theta) - model$p), rep(TRUE, model$p)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:(length(model$theta) - model$p),
                     pll = pppi, pman = pman, "ones", acut = acut, verbose = FALSE) #tape of the score function

  out <- smest(smoppi, 1:(length(model$theta) - model$p), model$sample, control = list(tol = 1E-10))

  expect_lt(out$value, smobj(smoppi, model$theta[1:(length(model$theta) - model$p)], model$sample))
  expect_lt(out$sqgradsize, sum(smobjgrad(smoppi, model$theta[1:(length(model$theta) - model$p)], model$sample)^2))

  cdabyppi:::expect_lt_v(abs(out$par - model$theta[1:(length(model$theta) - model$p)]) / out$SE, 3)
})

test_that("Fitting ppi via alr inc all beta gets close to true values", {
  set.seed(111)
  model <- sec2_3model(1000, maxden = 4)

  acut = 0.1 #not needed for Ralr per se, but code still expects it

  pman <- pmanifold("Ralr")
  pppi <- ptapell(rep(0.1, model$p - 1), model$theta, llname = "ppi", pman,
                  fixedtheta = rep(FALSE, length(model$theta)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:(length(model$theta)),
                     pll = pppi, pman = pman, "ones", acut = acut, verbose = FALSE) #tape of the score function

  out <- smest(smoppi, 1:length(model$theta), model$sample, control = list(tol = 1E-10))

  expect_lt(out$value, smobj(smoppi, model$theta, model$sample))
  expect_lt(out$sqgradsize, sum(smobjgrad(smoppi, model$theta, model$sample)^2))

  cdabyppi:::expect_lt_v(abs(out$par - model$theta) / out$SE, 3)
  # eestimating the -0.8 betas is poor
})
