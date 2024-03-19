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
})


test_that("tapell for vMF errors if NDEBUG not defined when theta isn't of the correct length", {
  # should get an assert error and R aborts
  maninfo <- manifoldtransform("sph")
  vmftape <- tapell(llname = "vMF",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = c(NA, NA),
                  tranobj = maninfo$tran)
})
