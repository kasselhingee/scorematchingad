test_that("tapell generates correct objects", {
maninfo <- manifoldtransform("sim", "sqrt", "sph")
ppitape <- tapell(llname = "ppi",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  tranobj = maninfo$tran) 
expect_true(R6::is.R6(ppitape))

# and verbose works too
expect_output({ppitape <- tapell(llname = "ppi",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  tranobj = maninfo$tran,
                  verbose = TRUE)},
              "pattern.*tape.*dynamic")
expect_true(R6::is.R6(ppitape))

expect_error({tapell(llname = "error",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  tranobj = maninfo$tran)}, "function")
})

test_that("tapell for ppi errors when theta isn't of the correct length", {
  maninfo <- manifoldtransform("sim", "sqrt", "sph")
  expect_condition(ppitape <- tapell(llname = "ppi",
                  ytape = c(0.2, 0.3, 0.3, 0.2),
                  usertheta = ppi_paramvec(p = 3),
                  tranobj = maninfo$tran), class = "Rcpp::exception", regexp = "length")

 # and that the taping was aborted too
 ppitape <- tapell(llname = "ppi",
                  ytape = c(0.2, 0.3, 0.3, 0.2),
                  usertheta = ppi_paramvec(p = 4),
                  tranobj = maninfo$tran)
})


test_that("tapell for vMF errors if NDEBUG not defined when theta isn't of the correct length", {
  skip("C++ assert errors not yet translated to R")
  # should get an assert error and R aborts
  maninfo <- manifoldtransform("sph")
  expect_condition(vmftape <- tapell(llname = "vMF",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = c(NA, NA),
                  tranobj = maninfo$tran))
})

test_that("llptr() internal function returns points like RcppXPtrUtils", {
  ppiXPtr <- getllptr("ppi")
  expect_type(ppiXPtr, "externalptr")


  #evaluate ppiXPtr (rather complicated at the moment)

psimplex <- manifoldtransform("sim", "identity", "sim")
lltape <- ptapell2(rep(1/3, 3), rep(-0.5, 3), llfXPtr = ppiXPtr, tran = psimplex$tran, fixedtheta = rep(FALSE, 3), verbose = FALSE)

expect_equal(3 * (-0.5 * log(1/3)), pForward0(lltape, rep(1/3, 3), rep(-0.5, 3)))
expect_equal(rep(-0.5 * 3, 3), pJacobian(lltape, rep(1/3, 3), rep(-0.5, 3)))


})
