test_that("Hess + Offset match gradient for a PPI Example", {
  mod <- ppi_egmodel(1)

  Ralr <- pmanifold("Ralr")
  ppitape <- tapell(llname = "ppi",
                  xtape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  pmanifoldtransform = Ralr)
  smotape <- tapesmo(lltape = ppitape,
                      pmanifoldtransform = Ralr,
                      divweight = "ones",
                      verbose = FALSE)

  testquadratictape(smotape)

  Hesstape <- pTapeHessian(smotape, attr(smotape, "xtape"), attr(smotape, "dyntape"))
  OffsetTape <- pTapeGradOffset(smotape, attr(smotape, "xtape"), attr(smotape, "dyntape"))

  #expect results to match for gradient
  x <- c(0.1, 0.1, 0.8)
  theta <- mod$theta
  grad1poly <- matrix(pForward0(Hesstape, 0*theta, x), ncol = length(theta)) %*% theta + 
    pForward0(OffsetTape, x, vector(mode = "double")) 
  grad1orig <- pJacobian(smotape, theta, x)
  expect_equal(grad1poly, grad1orig, ignore_attr = TRUE)

  # if W is symmetric then should equal the smo up to a constant
  # my expectation is that W is symmetric for ppi
  smo1poly <- 0.5 * theta %*% matrix(pForward0(Hesstape, 0*theta, x), ncol = length(theta)) %*% theta + 
    pForward0(OffsetTape, x, vector(mode = "double")) %*% theta

  x2 <- c(0.1, 0.4, 0.5)
  smo2poly <- 0.5 * theta %*% matrix(pForward0(Hesstape, 0*theta, x2), ncol = length(theta)) %*% theta + 
    pForward0(OffsetTape, x2, vector(mode = "double")) %*% theta

  smo1orig <- pForward0(smotape, theta, x)
  smo2orig <- pForward0(smotape, theta, x2)
  expect_equal(smo1poly - smo1orig, smo2poly - smo2orig, ignore_attr = TRUE)

}

