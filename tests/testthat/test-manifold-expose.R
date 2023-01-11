test_that("Manifold objects can be created", {
  mod <- Rcpp::Module("manif", PACKAGE="scorecompdir")
  manifold <- mod$manifold
  Ralr <- new(manifold, "Ralr")
})
