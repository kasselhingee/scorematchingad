
test_that("Printing CppAD graph of tape runs for dirichlet on simplex", {
  beta <- c(-0.8, -0.5, 1)
  pman <- pmanifold("simplex")
  pppi <- ptapell(rep(0.1, length(beta)), beta, llname = "dirichlet", pman,
                  fixedtheta = rep(FALSE, length(beta)), verbose = FALSE)
  expect_output(printgraph(pppi), regex = "mul.*log")
})

