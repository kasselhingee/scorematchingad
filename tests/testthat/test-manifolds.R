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
})

mod <- Rcpp::Module("manifolds", PACKAGE="scorecompdir")
obj <- mod$mantran_ad

test_that("Sphere manifold object matches analytic results", {
  sphere <- new(obj, "sphere")
  u <- c(0.1, 0.3, 0.6)
  z <- sphere$toM(u)
  expect_equal(sum(z^2), 1)
  expect_equal(sphere$fromM(z), u)
  Pmatz <- sphere$Pmatfun(z)
  expect_equal(t(z) %*% Pmatz %*% runif(3), 0, ignore_attr = TRUE)
  
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
  
  #check of first dPmatfun for first component
  z1 <- z[1]
  numgrad <- numericDeriv(quote(sphere$Pmatfun(c(z1, z[2:3]))), theta = "z1")
  expect_equal(drop(attr(numgrad, "gradient")), sphere$dPmatfun(z, 0))
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

test_that("Ralr manifold object matches analytic results",{
  Ralr <- new(obj, "Ralr")
  u <- c(0.1, 0.3, 0.6)
  z <- Ralr$toM(u)
  expect_length(z, 2)
  expect_equal(Ralr$fromM(z), u)
  expect_equal(Ralr$Pmatfun(z), diag(rep(1, 2)))
  expect_equal(Ralr$dPmatfun(z, 2), matrix(0, nrow = 2, ncol = 2))

  # check determinant
  integrand <- function(zmat){#each column is a measurement
    Jdets <- exp(apply(zmat, MARGIN = 2, Ralr$logdetJfromM))
    return(matrix(Jdets, nrow = 1))
  } 
  volume <- cubature::hcubature(
    f = integrand,
    lowerLimit = c(-1, -1) * 1E2,
    upperLimit = c(1, 1) * 1E2,
    vectorInterface = TRUE,
    fDim = 1)
  expect_equal(volume$integral, 0.5, tolerance = 1E-3) #0.5 seems to be the area of simplex under local coordinates (got 0.5 from integrating the simplex against local coordinates (e1, e2))
})

test_that("Snative manifold object matches analytic results", {
  sphere <- new(obj, "Snative")
  u <- runif(3)
  u <- u/sqrt(sum(u^2))
  expect_equal(u, sphere$toM(u))
  expect_equal(sphere$fromM(u), u)
  z <- u
  Pmatz <- sphere$Pmatfun(z)
  expect_equal(t(z) %*% Pmatz %*% (3*z), 0, ignore_attr = TRUE)
  expect_equal(t(z) %*% Pmatz %*% runif(3), 0, ignore_attr = TRUE)
  
  # determinant should be 1
  expect_equal(sphere$logdetJfromM(u), log(1))

  #check of first dPmatfun for first component
  z1 <- z[1]
  numgrad <- numericDeriv(quote(sphere$Pmatfun(c(z1, z[2:3]))), theta = "z1")
  expect_equal(drop(attr(numgrad, "gradient")), sphere$dPmatfun(z, 0))
})

test_that("Hclr manifold passes analytic tests",{
clr <- function(Y){
  logY <- log(Y)
  lgeommean <- rowSums(logY)/3
  return(logY - lgeommean)
}
clrinv <- function(Z){
  expZ <- exp(Z)
  return(expZ/rowSums(expZ))
}

ldetJfromM <- function(z){
  u <- as.vector(clrinv(matrix(z, nrow = 1)))
  sum(log(u)) + log(length(u))
}

  Hclr <- new(obj, "Hclr")
  u <- c(0.1, 0.3, 0.6)
  z <- Hclr$toM(u)
  expect_equal(z, clr(matrix(u, nrow = 1)), ignore_attr = TRUE)
  expect_equal(rep(1, 3) %*% z, 0, ignore_attr = TRUE)
  expect_equal(Hclr$fromM(z), u)

  # test determinant
  expect_equal(ldetJfromM(matrix(z, nrow = 1)), Hclr$logdetJfromM(z))

  # by integration
  

})



