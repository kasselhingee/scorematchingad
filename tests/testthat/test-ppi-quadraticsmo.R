## ppi() fitting tested with cppad_closed() tests

test_that("testquadratic passes on PPI model with sqrt transformation, minsq divergence weight, acut of 0.1", {
  tapes <- buildsmdtape(
     start = "sim",
     tran = "sqrt",
     end = "sph",
     ll = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3), 
     bdryw = "minsq",
     acut = 0.1,
     verbose = FALSE)
  ppismdtape <- tapes$smdtape


  # check only with parameter()
  expect_true(testquadratic(ppismdtape))



  xmat <- matrix(runif(length(ppi_paramvec(p = 3))*10, min = -0.9, max = 1),
                 ncol = length(ppi_paramvec(p = 3)))
  dynparammat <- matrix(runif(10*2, min=0, max = 0.5), ncol = 2)
  dynparammat <- cbind(dynparammat, 1-rowSums(dynparammat))

  expect_true(testquadratic(ppismdtape, xmat = xmat, dynparammat = dynparammat))
})

  # manual tests
test_that("manual tests on PPI model with sqrt transformation, minsq divergence weight, acut of 0.1", {
  tapes <- buildsmdtape(
     start = "sim",
     tran = "sqrt",
     end = "sph",
     ll = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3), 
     bdryw = "minsq",
     acut = 0.1,
     verbose = FALSE)
  ppismdtape <- tapes$smdtape

  ppismdtapeJ <- tapeJacobian(ppismdtape)
  ppismdtapeH <- tapeHessian(ppismdtape)
  ppismdtapeH2 <- tapeJacobian(ppismdtapeJ)

  expect_equal(
  ppismdtapeH$eval(ppismdtape$xtape, c(0.1, 0.1, 0.8)),
  ppismdtapeH2$eval(ppismdtape$xtape, c(0.1, 0.1, 0.8)))

  #test that the value of Hessian *does* depend on the measurement vector  
  Hnearedge <- ppismdtapeH2$eval(ppismdtape$xtape, c(0.1, 0.1, 0.8))
  Hnearcentre <- ppismdtapeH2$eval(ppismdtape$xtape, c(0.3, 0.3, 0.4))
  expect_false(isTRUE(all.equal(Hnearedge, Hnearcentre)))

  # the next results are false for a reason unknown to me because the tape seems to be doing the right thing
  expect_equal(sapply(1:ppismdtapeH$range, function(i){ppismdtapeH$parameter(i-1)}), rep(FALSE, length(ppi_paramvec(p = 3))^2))
  expect_equal(sapply(1:ppismdtapeH2$range, function(i){ppismdtapeH2$parameter(i-1)}), rep(TRUE, length(ppi_paramvec(p = 3))^2))
  
  hessgrad <- ppismdtapeH$Jac(
                           ppi_paramvec(p = 3, AL=1, bL=1, beta=c(-0.1,-0.1,0.5)),
                           c(0.1, 0.1, 0.8))
  expect_equal(hessgrad, rep(0, length(ppi_paramvec(p = 3))^3))
})

test_that("ppi ll tape is fails the quadratic test", {
  sqrtman <- manifoldtransform("sim", "sqrt", "sph")
  ppitape <- tapell(ll = "ppi",
                    ytape = c(0.2, 0.3, 0.5),
                    usertheta = ppi_paramvec(p = 3), 
                    tranobj = sqrtman$tran)

  # check only with pParameter()
  expect_false(testquadratic(ppitape))
  expect_message(testquadratic(ppitape, verbose = TRUE), "non-constant")

  # check some values

  dynparammat <- matrix(runif(length(ppi_paramvec(p = 3))*10, min = -0.9, max = 1),
                 ncol = length(ppi_paramvec(p = 3)))
  xmat <- matrix(runif(10*2, min=0, max = 0.5), ncol = 2)
  xmat <- cbind(xmat, 1-rowSums(xmat))

  expect_false(testquadratic(ppitape, xmat = xmat, dynparammat = dynparammat))
  testquadratic(ppitape, xmat = xmat, dynparammat = dynparammat, verbose = TRUE) |> 
    expect_message("non-constant") |>
    expect_message("non-zero")
})

