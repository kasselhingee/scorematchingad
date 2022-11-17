test_that("Solution without boundary considerations for PPI has zero gradient and matches numerical minimum", {
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
