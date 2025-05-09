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

test_that("tape_bdryw_inbuilt() works", {
  tape <- tape_bdryw_inbuilt("minsq", rep(0.2, 5), 0.1)
  expect_equal(tape$xtape, rep(0.2, 5))

  expect_equal(tape$forward(0, rep(0.2, 5)), 0.1^2)
  expect_equal(tape$forward(0, c(0.01, 0.19, 0.2, 0.2, 0.2)), 0.01^2)

  expect_equal(tape$Jacobian(rep(0.2, 5)), rep(0, 5))
  expect_equal(tape$Jacobian(c(0.01, 0.19, 0.2, 0.2, 0.2)), c(2*0.01, rep(0, 4)))
})
