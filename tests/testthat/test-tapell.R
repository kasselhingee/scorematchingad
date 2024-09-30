test_that("tapell generates correct objects", {
maninfo <- manifoldtransform("sim", "sqrt", "sph")
ppitape <- tapell(ll = "ppi",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  tranobj = maninfo$tran)
expect_s4_class(ppitape, "Rcpp_ADFun")


expect_error({tapell(ll = "error",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  tranobj = maninfo$tran)}, "function")
})


test_that("tapell for ppi errors when theta isn't of the correct length", {
  maninfo <- manifoldtransform("sim", "sqrt", "sph")
  expect_condition(ppitape <- tapell(ll = "ppi",
                  ytape = c(0.2, 0.3, 0.3, 0.2),
                  usertheta = ppi_paramvec(p = 3),
                  tranobj = maninfo$tran), class = "Rcpp::exception", regexp = "length")

 # and that the taping was aborted too
 ppitape <- tapell(ll = "ppi",
                  ytape = c(0.2, 0.3, 0.3, 0.2),
                  usertheta = ppi_paramvec(p = 4),
                  tranobj = maninfo$tran)
})

