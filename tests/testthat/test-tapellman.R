test_that("tapellman() has correct range and domain length, and gives correct values", {
maninfo <- manifoldtransform("sim", "identity", "sim")
lltape <- tapell(ll = "dirichlet",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = rep(NA, 3),
                  tranobj = maninfo$tran)

maninfo <- manifoldtransform("sim", "sqrt", "sph")
llmantape <- ptapellman(c(0.2, 0.3, 0.5), rep(1, 3), lltape$ptr, maninfo$tran, verbose = FALSE)

expect_equal(pForward0(llmantape, rep(1/3, 3), rep(-0.5, 3)), 2.0794416)
expect_equal(pJacobian(llmantape, c(0.1, 0.1, 0.8), rep(-0.1, 3)), c(8,8,1))
})

