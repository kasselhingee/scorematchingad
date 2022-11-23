test_that("testquadratictape passes on PPI model with sqrt transformation, minsq divergence weight, acut of 0.1", {
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

  # check only with pParameter()
  expect_true(testquadratictape(ppismotape))



  xmat <- matrix(runif(length(ppi_paramvec(p = 3))*10, min = -0.9, max = 1),
                 ncol = length(ppi_paramvec(p = 3)))
  dynparammat <- matrix(runif(10*2, min=0, max = 0.5), ncol = 2)
  dynparammat <- cbind(dynparammat, 1-rowSums(dynparammat))

  expect_true(testquadratictape(ppismotape, xmat = xmat, dynparammat = dynparammat))
})

  # manual tests
test_that("manual tests on PPI model with sqrt transformation, minsq divergence weight, acut of 0.1", {
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

  #test that the value of Hessian *does* depend on the measurement vector  
  Hnearedge <- pForward0(ppismotapeH2, attr(ppismotape, "xtape"), c(0.1, 0.1, 0.8))
  Hnearcentre <- pForward0(ppismotapeH2, attr(ppismotape, "xtape"), c(0.3, 0.3, 0.4))
  expect_false(isTRUE(all.equal(Hnearedge, Hnearcentre)))

  # the next results are false for a reason unknown to me because the tape seems to be doing the right thing
  expect_equal(pParameter(ppismotapeH), rep(FALSE, length(ppi_paramvec(p = 3))^2))
  expect_equal(pParameter(ppismotapeH2), rep(TRUE, length(ppi_paramvec(p = 3))^2))
  
  hessgrad <- pJacobian(ppismotapeH,
                           ppi_paramvec(p = 3, AL=1, bL=1, beta=c(-0.1,-0.1,0.5)),
                           c(0.1, 0.1, 0.8))
  expect_equal(hessgrad, rep(0, length(ppi_paramvec(p = 3))^3))
})

test_that("ppi ll tape is fails the quadratic test", {
  sqrtman <- pmanifold("sphere")
  ppitape <- tapell(llname = "ppi",
                    xtape = c(0.2, 0.3, 0.5),
                    usertheta = ppi_paramvec(p = 3), 
                    pmanifoldtransform = sqrtman)

  # check only with pParameter()
  expect_false(testquadratictape(ppitape))
  expect_message(testquadratictape(ppitape, verbose = TRUE))

  # check some values

  dynparammat <- matrix(runif(length(ppi_paramvec(p = 3))*10, min = -0.9, max = 1),
                 ncol = length(ppi_paramvec(p = 3)))
  xmat <- matrix(runif(10*2, min=0, max = 0.5), ncol = 2)
  xmat <- cbind(xmat, 1-rowSums(xmat))

  expect_false(testquadratictape(ppitape, xmat = xmat, dynparammat = dynparammat))
  expect_message(testquadratictape(ppitape, xmat = xmat, dynparammat = dynparammat, verbose = TRUE))
})

