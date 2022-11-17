test_that("Hess + Offset match gradient for a PPI Example", {
  mod <- ppi_egmodel(100)
  Y <- mod$sample

  Ralr <- pmanifold("Ralr")
  ppitape <- tapell(llname = "ppi",
                  xtape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  pmanifoldtransform = Ralr)
  smotape <- tapesmo(lltape = ppitape,
                      pmanifoldtransform = Ralr,
                      divweight = "ones",
                      verbose = FALSE)

  values <- quadratictape_parts(smotape, Y)

  # expect results to match for gradient
  theta <- mod$theta #theta could be anything though
  gradorig <- t(apply(Y, MARGIN = 1, function(x){pJacobian(smotape, theta, x)}))

  gradpoly <- lapply(1:nrow(values$offset), function(i){
    drop(matrix(values$Hessian[i, ], ncol = length(theta)) %*% theta + 
      t(values$offset[i, , drop = FALSE]))
  })
  gradpoly <- do.call(rbind, gradpoly)
  expect_equal(gradorig, gradpoly)

  # if W is symmetric then should equal the smo up to a constant
  smoorig <- apply(Y, MARGIN = 1, function(x){pForward0(smotape, theta, x)})

  # my expectation is that W is symmetric for ppi
  smopoly <- lapply(1:nrow(values$offset), function(i){
    drop(0.5 * theta %*% matrix(values$Hessian[i, ], ncol = length(theta)) %*% theta + 
      values$offset[i, , drop = FALSE] %*% theta)
  })
  smopoly <- unlist(smopoly)

  expect_true(all(smoorig - smopoly == smoorig[[1]] - smopoly[[1]])) #seems I'm wrong!? - because there is a free parameter in the way I've set up ppi
}

test_that("quadratictape_parts_approx is close to quadratic_parts for simplex interior points", {
  mod <- ppi_egmodel(100)
  Y <- mod$sample
  Ycen <- simplex_boundaryshift(Y)

  Ralr <- pmanifold("Ralr")
  ppitape <- tapell(llname = "ppi",
                  xtape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  pmanifoldtransform = Ralr)
  smotape <- tapesmo(lltape = ppitape,
                      pmanifoldtransform = Ralr,
                      divweight = "ones",
                      verbose = FALSE)
  
  valuesexact <- quadratictape_parts(smotape, Y)
  valuesapprox <- quadratictape_parts_approx(smotape, Y, centres = Ycen, order = 10)

})
