test_that("Fitting ppi via alr transform gets close to true values", {
  set.seed(111)
  model <- sec2_3model(1000, maxden = 4)

  acut = 0.1 #not needed for Ralr per se, but code still expects it

  pman <- pmanifold("Ralr")
  pppi <- ptapell(rep(0.1, model$p - 1), model$theta, llname = "ppi", pman, fixedtheta = c(rep(FALSE, length(model$theta) - 1), TRUE), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, model$p), 1:(length(model$theta) - 1), pll = pppi, pman = pman, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  out <- smest(smoppi, 1:(length(model$theta) - 1), model$sample, control = list(tol = 1E-20))

  cdabyppi:::expect_lt_v(abs(out$par - model$theta[-length(model$theta)]) / out$SE, 3)
})
