test_that("tape_smi() return value is correct", {
  # set.seed(123)
  m <- rppi_egmodel(2)
  acut <- 0.1
  
  tapes <- tape_smi(manifold = "sph",
                    uld = "ppi",
                    transform = "sqrt",
                    xtape = c(0.1,0.1,0.1), 
                    fixedparams = rep(NA, length(m$theta)),
                    bdryw =  tape_bdryw_inbuilt("minsq", c(0.1,0.1,0.1), acut = acut))
  expect_equal(tapes$info$transform, "sqrt")
  expect_equal(tapes$info$manifold, "sph")
  expect_s4_class(tapes$info$bdryw, "Rcpp_ADFun")
  expect_s4_class(tapes$smi, "Rcpp_ADFun")
  expect_s4_class(tapes$uld_reembed, "Rcpp_ADFun")
  
  # check properties of returned tapes
  expect_equal(tapes$smi$domain, length(m$theta))
  expect_equal(tapes$smi$range, 1)
  expect_equal(tapes$smi$size_dyn_ind, 3)
  expect_equal(tapes$uld_reembed$domain, 3)
  expect_equal(tapes$uld_reembed$range, 1)
  expect_equal(tapes$uld_reembed$size_dyn_ind, length(m$theta))
  expect_equal(tapes$uld_reembed$xtape, sqrt(c(0.1, 0.1, 0.1)))
})






