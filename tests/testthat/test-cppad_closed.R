test_that("Solution without boundary considerations for PPI has zero gradient and matches numerical minimum", {
  set.seed(13411)
  mod <- ppi_egmodel(100)
  Y <- mod$sample

  Ralr <- pmanifold("Ralr")
  ppitape <- tapell(llname = "ppi",
                  xtape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3, betap = tail(mod$beta, 1)), 
                  pmanifoldtransform = Ralr)
  smotape <- tapesmo(lltape = ppitape,
                      pmanifoldtransform = Ralr,
                      divweight = "ones",
                      verbose = FALSE)

  estobj <- cppad_closed(smotape, Y)

  grads <- t(apply(Y, MARGIN = 1, function(x) pJacobian(smotape, estobj$est, x))) 
  totalgrad <- colSums(grads)
  expect_lt(sum(totalgrad^2), 1E-20)

  numericalmin <- ppi(Y, paramvec = ppi_paramvec(p = 3, betap = tail(mod$beta, 1)), trans = "alr", method = "cppad")
  expect_equal(numericalmin$est$paramvec, c(estobj$est, tail(mod$beta, 1)), ignore_attr = TRUE)
})

test_that("Closed-from solution with boundary points matches hard-coded version", {
  mnongamma <- ppi_egmodel(1)
  theta <- ppi_paramvec(beta = c(-0.95, -0.9, 0.5), AL = mnongamma$AL, bL = 0)
  set.seed(1234)
  Ycts <- rppi(1000, paramvec = theta)
  dsample <- round(Ycts * 100)/ 100
  dsample[, 3] <- 1 - rowSums(dsample[, 1:2])

  isbdry <- simplex_isboundary(dsample)
  Yapproxcentres <- dsample
  Yapproxcentres[!isbdry, ] <- NA 
  Yapproxcentres[isbdry, ] <- simplex_boundaryshift(dsample[isbdry, , drop = FALSE])

  Ralr <- pmanifold("Ralr")
  ppitape <- tapell(llname = "ppi",
                  xtape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3, bL = 0, betap = tail(theta, 1)), 
                  pmanifoldtransform = Ralr)
  smotape <- tapesmo(lltape = ppitape,
                      pmanifoldtransform = Ralr,
                      divweight = "ones",
                      verbose = FALSE)

  estobj <- cppad_closed(smotape, Y = dsample, Yapproxcentres, approxorder = 10)

  est_hardcode <- ppi(dsample, paramvec = ppi_paramvec(p = 3, bL = 0, betap = tail(theta, 1)),
      trans = "alr", method = "direct")
  expect_equal(est_hardcode$est$paramvec, t_fu2t(estobj$est, attr(ppitape, "usertheta")))

})
