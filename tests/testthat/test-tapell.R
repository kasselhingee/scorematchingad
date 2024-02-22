test_that("tapell generates correct objects", {
maninfo <- manifoldtransform("sim", "sqrt", "sph")
ppitape <- tapell(llname = "ppi",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  tran = maninfo$tran) 
expect_true(R6::is.R6(ppitape))

# and verbose works too
expect_output({ppitape <- tapell(llname = "ppi",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  tran = maninfo$tran,
                  verbose = TRUE)},
              "pattern.*tape.*dynamic")
expect_true(R6::is.R6(ppitape))

expect_error({tapell(llname = "error",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  tran = maninfo$tran)}, "function")
})
