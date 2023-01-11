test_that("Manifold objects can be created", {
  mod <- Rcpp::Module("manifolds", PACKAGE="scorecompdir")
  mantranobj <- mod$mantranobj
  Ralr <- new(mantranobj, "Ralr")
  Ralr$toM(c(0.1, 0.3, 0.6))

})
