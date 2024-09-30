test_that("reembed() is correct for sqrt transform of simplex metric", {
  utape = rep(0.2, 5)
  dyntape = rep(-0.1, 5)
  tape <- tape_uld_inbuilt("dirichlet", utape, dyntape)


  maninfo <- manifoldtransform("sim", "sqrt", "sph")
  dirichwrtsph <- reembed(tape, maninfo$tran)
  expect_equal(dirichwrtsph$xtape, sqrt(utape))
  expect_equal(dirichwrtsph$dyntape, dyntape)

  newu <- c(0.1, 0.1, 0.2, 0.2, 0.4)
  newbeta <- c(-0.9, -0.5, -0.1, -0.1, -0.1)
  expect_equal(dirichwrtsph$Jac(sqrt(newu), newbeta), (1+2*newbeta)/sqrt(newu))
})


