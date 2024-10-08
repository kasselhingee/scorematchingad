test_that("ADFun object is returned by tapell, and its values can be accessed", {
maninfo <- manifoldtransform("sim", "sqrt", "sph")
ppitape <- tapell(ll = "ppi",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  tranobj = maninfo$tran)

expect_type(ppitape$xtape, "double") #vectors are neither S3 objects nor S4 objects
expect_type(ppitape$dyntape, "double")
expect_type(ppitape$name, "character")
})

test_that("print.Rcpp_ADFun() works", {
  utape = rep(0.2, 5)
  dyntape = rep(-0.1, 5)
  tape <- tape_uld_inbuilt("dirichlet", utape, dyntape)
  expect_output(print(tape), "dirichlet.*5")
})
