test_that("Hess + Offset match gradient for a PPI Example", {
  mod <- ppi_egmodel(100)
  Y <- mod$sample

  tapes <- buildsmotape(
     manifoldname = "Ralr",
     llname = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3, betap=0.5), 
     verbose = FALSE)
  smotape <- tapes$smotape

  values <- quadratictape_parts(smotape, Y)

  # expect results to match for gradient
  theta <- head(mod$theta, 7) #theta could be anything though
  gradorig <- t(apply(Y, MARGIN = 1, function(x){pJacobian(smotape$ptr, theta, x)}))

  gradpoly <- lapply(1:nrow(values$offset), function(i){
    drop(matrix(values$Hessian[i, ], ncol = length(theta)) %*% theta + 
      t(values$offset[i, , drop = FALSE]))
  })
  gradpoly <- do.call(rbind, gradpoly)
  expect_equal(gradorig, gradpoly)

  # if W is symmetric then should equal the smo up to a constant wrt theta
  # and my expectation is that W is symmetric for ppi
  smoorig <- evaltape(smotape, theta, Y)

  smopoly <- lapply(1:nrow(values$offset), function(i){
    drop(0.5 * theta %*% matrix(values$Hessian[i, ], ncol = length(theta)) %*% theta + 
      values$offset[i, , drop = FALSE] %*% theta)
  })
  smopoly <- unlist(smopoly)
  constant <- smoorig-smopoly
  
  #test constant by trying another theta
  theta2 <- theta+1
  smoorig2 <- apply(Y, MARGIN = 1, function(x){pForward0(smotape, theta2, x)})
  smopoly2 <- lapply(1:nrow(values$offset), function(i){
    drop(0.5 * theta2 %*% matrix(values$Hessian[i, ], ncol = length(theta2)) %*% theta2 +
      values$offset[i, , drop = FALSE] %*% theta2)
  })
  smopoly2 <- unlist(smopoly2)
  expect_equal(smoorig2-smopoly2, constant)
})

test_that("quadratictape_parts with approx centres is close to quadratic_parts for simplex interior points", {
  mod <- ppi_egmodel(100)
  Y <- mod$sample
  Ycen <- simplex_boundaryshift(Y)

  tapes <- buildsmotape(
     manifoldname = "Ralr",
     llname = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3), 
     verbose = FALSE)
  smotape <- tapes$smotape
  
  valuesexact <- quadratictape_parts(smotape, Y)
  valuesapprox <- quadratictape_parts(smotape, Y, tcentres = Ycen, approxorder = 1)
  expect_equal(valuesexact, valuesapprox, tolerance = 1E-5)
  # but still an approximation
  expect_gt(max(abs(valuesexact$offset - valuesapprox$offset)), 1E-10)

  #expect higher order approximation to be closer
  valuesapprox2 <- quadratictape_parts(smotape, Y, tcentres = Ycen, approxorder = 10)
  expect_lt(sum((valuesexact$offset - valuesapprox2$offset)^2), 
            sum((valuesexact$offset - valuesapprox$offset)^2))
})


