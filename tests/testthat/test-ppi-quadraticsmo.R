## ppi() fitting tested with cppad_closed() tests

test_that("testquadratictape passes on PPI model with sqrt transformation, minsq divergence weight, acut of 0.1", {
  tapes <- buildsmotape(
     manifoldname = "sphere",
     llname = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3), 
     weightname = "minsq",
     acut = 0.1,
     verbose = FALSE)
  ppismotape <- tapes$smotape


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
  tapes <- buildsmotape(
     manifoldname = "sphere",
     llname = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3), 
     weightname = "minsq",
     acut = 0.1,
     verbose = FALSE)
  ppismotape <- tapes$smotape

  ppismotapeJ <- tapeJacobian(ppismotape)
  ppismotapeH <- tapeHessian(ppismotape)
  ppismotapeH2 <- tapeJacobian(ppismotapeJ)

  expect_equal(
  pForward0(ppismotapeH$ptr, ppismotape$xtape, c(0.1, 0.1, 0.8)),
  pForward0(ppismotapeH2$ptr, ppismotape$xtape, c(0.1, 0.1, 0.8)))

  #test that the value of Hessian *does* depend on the measurement vector  
  Hnearedge <- pForward0(ppismotapeH2$ptr, ppismotape$xtape, c(0.1, 0.1, 0.8))
  Hnearcentre <- pForward0(ppismotapeH2$ptr, ppismotape$xtape, c(0.3, 0.3, 0.4))
  expect_false(isTRUE(all.equal(Hnearedge, Hnearcentre)))

  # the next results are false for a reason unknown to me because the tape seems to be doing the right thing
  expect_equal(pParameter(ppismotapeH$ptr), rep(FALSE, length(ppi_paramvec(p = 3))^2))
  expect_equal(pParameter(ppismotapeH2$ptr), rep(TRUE, length(ppi_paramvec(p = 3))^2))
  
  hessgrad <- pJacobian(ppismotapeH$ptr,
                           ppi_paramvec(p = 3, AL=1, bL=1, beta=c(-0.1,-0.1,0.5)),
                           c(0.1, 0.1, 0.8))
  expect_equal(hessgrad, rep(0, length(ppi_paramvec(p = 3))^3))
})

test_that("ppi ll tape is fails the quadratic test", {
  sqrtman <- manifoldtransform("sphere")
  ppitape <- tapell(llname = "ppi",
                    ytape = c(0.2, 0.3, 0.5),
                    usertheta = ppi_paramvec(p = 3), 
                    pmanifoldtransform = sqrtman)

  # check only with pParameter()
  expect_false(testquadratictape(ppitape$ptr))
  expect_message(testquadratictape(ppitape$ptr, verbose = TRUE), "non-constant")

  # check some values

  dynparammat <- matrix(runif(length(ppi_paramvec(p = 3))*10, min = -0.9, max = 1),
                 ncol = length(ppi_paramvec(p = 3)))
  xmat <- matrix(runif(10*2, min=0, max = 0.5), ncol = 2)
  xmat <- cbind(xmat, 1-rowSums(xmat))

  expect_false(testquadratictape(ppitape, xmat = xmat, dynparammat = dynparammat))
  testquadratictape(ppitape, xmat = xmat, dynparammat = dynparammat, verbose = TRUE) |> 
    expect_message("non-constant") |>
    expect_message("non-zero")
})

