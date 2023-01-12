test_that("Manifold objects can be created, and member functions run", {
  mod <- Rcpp::Module("manifolds", PACKAGE="scorecompdir")
  obj <- mod$mantran_ad
  Ralr <- new(obj, "Ralr")
  expect_s4_class(Ralr, "Rcpp_mantran_ad")
  u <- c(0.1, 0.3, 0.6)
  z <- Ralr$toM(u)
  expect_equal(Ralr$fromM(z), u)
  expect_true(is.finite(Ralr$logdetJfromM(z)))
  expect_equal(dim(Ralr$Pmatfun(z)), c(2,2))
  expect_type(Ralr$Pmatfun(z), "double")
  expect_equal(dim(Ralr$dPmatfun(z, 1)), c(2,2))
  expect_type(Ralr$XPtr(), "externalptr")
})


