test_that("Fitting ppi via clr transform with fixed beta gets close to true values", {
  skip_on_cran()
  set.seed(1234)
  model <- ppi_egmodel(1000, maxden = 4)

  out <- ppi(Y = model$sample,
             paramvec = ppi_paramvec(beta = model$beta),
             trans = "clr")

  expect_absdiff_lte_v(out$est$paramvec, model$theta, out$SE$paramvec * 3)
})


test_that("Hess + Offset match gradient for Hclr", {
  mod <- ppi_egmodel(100)
  Y <- mod$sample

  Hclr <- manifoldtransform("Hclr")
  ppitape <- tapell(llname = "ppi",
                  xtape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  pmanifoldtransform = Hclr)

  tape_eval(ppitape, Y[2, ], mod$theta)
  tape_eval(ppitape, c(0.120937, 0.1082838, 0.7707792), mod$theta) #small change results in finite numerical evaluation - I think the way the Jacobian is differentiated is creating the problems.

  smotape <- tapesmo(lltape = ppitape,
                      pmanifoldtransform = Hclr,
                      divweight = "ones",
                      verbose = FALSE)


  values <- quadratictape_parts(smotape, Y)

  # expect results to match for gradient
  gradorig <- t(apply(Y, MARGIN = 1, function(x){pJacobian(smotape, theta, x)}))

  gradpoly <- lapply(1:nrow(values$offset), function(i){
    drop(matrix(values$Hessian[i, ], ncol = length(theta)) %*% theta + 
      t(values$offset[i, , drop = FALSE]))
  })
  gradpoly <- do.call(rbind, gradpoly)
  expect_equal(gradorig, gradpoly)

  # if W is symmetric then should equal the smo up to a constant wrt theta
  # and my expectation is that W is symmetric for ppi
  smoorig <- apply(Y, MARGIN = 1, function(x){pForward0(smotape, theta, x)})

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

