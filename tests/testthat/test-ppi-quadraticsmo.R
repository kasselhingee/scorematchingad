test_that("PPI model with sqrt transformation, minsq divergence weight, acut of 0.1 is quadratic", {
  sqrtman <- pmanifold("sphere")
  ppitape <- tapell(llname = "ppi",
                    xtape = c(0.2, 0.3, 0.5),
                    usertheta = ppi_paramvec(p = 3), 
                    pmanifoldtransform = sqrtman)
  ppismotape <- tapesmo(lltape = ppitape,
                        pmanifoldtransform = sqrtman,
                        divweight = "minsq",
                        acut = 0.1,
                        verbose = FALSE)
  ppismotapeJ <- pTapeJacobian(ppismotape, attr(ppismotape, "xtape"), attr(ppismotape, "dyntape"))
  ppismotapeH <- pTapeHessian(ppismotape, attr(ppismotape, "xtape"), attr(ppismotape, "dyntape"))
  ppismotapeH2 <- pTapeJacobian(ppismotapeJ, attr(ppismotape, "xtape"), attr(ppismotape, "dyntape"))

  expect_equal(
  pForward0(ppismotapeH, attr(ppismotape, "xtape"), c(0.1, 0.1, 0.8)),
  pForward0(ppismotapeH2, attr(ppismotape, "xtape"), c(0.1, 0.1, 0.8)))
  
  expect_equal(pParameter(ppismotapeH, c(0.1, 0.1, 0.8)), rep(TRUE, length(ppi_paramvec(p = 3))^2))
  expect_equal(pParameter(ppismotapeH2, c(0.1, 0.1, 0.8)), rep(TRUE, length(ppi_paramvec(p = 3))^2))
  
  hessgrad <- scorecompdir:::pJacobian(ppismotapeH,
                           ppi_paramvec(p = 3, AL=1, bL=1, beta=c(-0.1,-0.1,0.5)),
                           c(0.1, 0.1, 0.8))
  expect_equal(hessgrad, rep(0, length(ppi_paramvec(p = 3))^3))
})


