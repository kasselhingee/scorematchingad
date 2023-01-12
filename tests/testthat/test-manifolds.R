skip_on_cran()
mod <- Rcpp::Module("manifolds", PACKAGE="scorecompdir")
obj <- mod$mantran_ad

test_that("Sphere manifold passes lightweight standard tests",{
  pman <- manifoldtransform("sphere")
  u <- c(0.1, 0.2, 1 - 0.3)
  out <- testmanifold(pman, u)
  expect_equal(out, 0)
})

test_that("Sphere manifold object matches analytic results", {
  sphere <- new(obj, "sphere")
  u <- c(0.1, 0.3, 0.6)
  z <- sphere$toM(u)
  expect_equal(sum(z^2), 1)
  expect_equal(sphere$fromM(z), u)
  Pmatz <- sphere$Pmatfun(z)
  expect_equal(t(z) %*% Pmatz %*% runif(3), 0)
  
  #check determinant using integration over a unitcube
  integrand <- function(zmat){#each column is a measurement
    xmat <- apply(zmat, MARGIN = 2, sphere$fromM)
    inunitcube <- (colSums(xmat < 0) == 0) * (colSums(xmat > 1) == 0)
    Jdets <- exp(apply(zmat, MARGIN = 2, sphere$logdetJfromM))
    return(matrix(Jdets * inunitcube, nrow = 1))
  } 
  volumeviaM <- cubature::hcubature(
    f = integrand,
    lowerLimit = c(0, 0),
    upperLimit = c(1, 1),
    vectorInterface = TRUE,
    fDim = 1)
  expect_equal(volumeviaM$integral, 1, tolerance = 1E-5)
})

test_that("Simplex manifold object matches analytic results", {
  simplex <- new(obj, "simplex")
  u <- c(0.1, 0.3, 0.6)
  expect_equal(u, simplex$toM(u))
  expect_equal(u, simplex$fromM(u))
  expect_equal(simplex$logdetJfromM(u), log(1))
  expect_equal(simplex$Pmatfun(u), diag(rep(1, 3)) - rep(1, 3) %*% t(rep(1, 3))/3)
  # check projection
  set.seed(1342)
  x <- runif(3)
  u2 <- simplex$Pmatfun(u) %*% x
  expect_equal(t(u2) %*% rep(1, 3), 0, ignore_attr = TRUE)

  expect_equal(simplex$dPmatfun(u, 2), matrix(0, nrow = 3, ncol = 3))
})

test_that("Ralr manifold passes lightweight standard tests",{
  pman <- manifoldtransform("Ralr")
  u <- c(0.1, 0.2, 1 - 0.3)
  out <- testmanifold(pman, u)
  expect_equal(out, 0)
})

test_that("Snative manifold object passes lightweight standard tests", {
  pman <- manifoldtransform("Snative")
  u <- c(0.1, 0.2)
  u <- u / sqrt(sum(u^2))
  out <- testmanifold(pman, u)
  expect_equal(out, 0)
})

