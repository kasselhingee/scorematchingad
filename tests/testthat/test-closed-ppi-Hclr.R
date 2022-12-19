test_that("Fitting ppi via clr transform with fixed beta gets close to true values and other estimators", {
  skip_on_cran()
  set.seed(1234)
  model <- ppi_egmodel(1000, maxden = 4)

  out <- ppi(Y = model$sample,
             paramvec = ppi_paramvec(beta = model$beta),
             trans = "clr")
  expect_absdiff_lte_v(out$est$paramvec, model$theta, out$SE$paramvec * 3)
  # but the estimate is wrong by 10 orders of magnitude


  out2 <- ppi(Y = model$sample,
             paramvec = ppi_paramvec(beta = model$beta),
             trans = "sqrt",
             acut = 0.01, divweight = "minsq")

  expect_absdiff_lte_v(out$est$paramvec, out2$est$paramvec, out2$SE$paramvec * 3)
})
test_that("Fitting ppi via clr transform with unfixed beta gets close to true values", {
  skip_on_cran()
  set.seed(1234)
  model <- ppi_egmodel(1000, maxden = 4)

  out <- ppi(Y = model$sample,
             trans = "clr")

  expect_absdiff_lte_v(out$est$paramvec, model$theta, out$SE$paramvec * 3)
})

test_that("tapefromM evaluates, derivative and Jacobian match R-computed versions", {
  set.seed(1245)
  mod <- ppi_egmodel(100)
  Y <- mod$sample
  Y <- rbind(Y, rep(1/3, 3))

  Hclr <- manifoldtransform("Hclr")

  # check evaluation
  clrY <- log(Y) - rowSums(log(Y))/3
  tapefromM <- ptapefromM(clrY[1, ], Hclr)
  vals <- tape_eval(tapefromM, clrY, matrix(nrow = nrow(Y), ncol = 0))
  expect_equal(vals, Y)

  #check nondegeneracy fix
  vals2 <- tape_eval(tapefromM, clrY+3.1, matrix(nrow = nrow(Y), ncol = 0))
  expect_equal(vals2-3.1, Y)

  # check Jacobian
  fromMJ <- t(apply(clrY, MARGIN = 1, function(x){pJacobian(tapefromM, x, matrix(nrow = 0, ncol = 0))}))

  # numerical Jacobian
  fromMJnum <- t(apply(clrY, MARGIN = 1, function(x){
    b <- x
    J <- numericDeriv(quote(pForward0(tapefromM, b , matrix(nrow = 0, ncol = 0))), "b")
    return(as.vector(attr(J, "gradient")))
}))
  expect_equal(fromMJ, fromMJnum, tolerance = 1E-4)

  detJ <- apply(fromMJ, MARGIN = 1, function(J){
    Jmat <- matrix(J, nrow = sqrt(length(J)), ncol = sqrt(length(J)))
    return(det(Jmat))
  })

  detJnum <- apply(fromMJnum, MARGIN = 1, function(J){
    Jmat <- matrix(J, nrow = sqrt(length(J)), ncol = sqrt(length(J)))
    return(det(Jmat))
  })
  expect_equal(detJ, detJnum, tolerance = 1E-4)
})

test_that("Hess + Offset match gradient for Hclr", {
  set.seed(1245)
  mod <- ppi_egmodel(100)
  Y <- mod$sample
  Y <- rbind(Y, rep(1/3, 3))

  tapes <- buildsmotape("Hclr", "ppi", utape = c(0.2, 0.3, 0.5), 
                        usertheta = ppi_paramvec(p = 3))

  values <- quadratictape_parts(tapes$smotape, Y)

  # expect results to match for gradient
  gradorig <- t(apply(Y, MARGIN = 1, function(x){pJacobian(tapes$smotape, mod$theta, x)}))

  gradpoly <- lapply(1:nrow(values$offset), function(i){
    drop(matrix(values$Hessian[i, ], ncol = length(mod$theta)) %*% mod$theta + 
      t(values$offset[i, , drop = FALSE]))
  })
  gradpoly <- do.call(rbind, gradpoly)
  expect_equal(gradorig, gradpoly)
})

test_that("W is symmetric for ppi with clr, beta fixed", {
  # if W is symmetric then should equal the smo up to a constant wrt theta
  # and my expectation is that W is symmetric for ppi with Hclr
  set.seed(1245)
  mod <- ppi_egmodel(100)
  Y <- mod$sample
  Y <- rbind(Y, rep(1/3, 3))

  usertheta = ppi_paramvec(p = 3, beta = mod$beta)
  ftheta <- t_ut2f(usertheta, mod$theta)
  tapes <- buildsmotape("Hclr", "ppi", utape = c(0.2, 0.3, 0.5), 
                        usertheta = usertheta)

  values <- quadratictape_parts(tapes$smotape, Y)

  smoorig <- apply(Y, MARGIN = 1, function(x){pForward0(tapes$smotape, ftheta, x)})

  smopoly <- lapply(1:nrow(values$offset), function(i){
    drop(0.5 * ftheta %*% matrix(values$Hessian[i, ], ncol = length(ftheta)) %*% ftheta + 
      values$offset[i, , drop = FALSE] %*% ftheta)
  })
  smopoly <- unlist(smopoly)
  constant <- smoorig-smopoly
  
  #test constant by trying another theta
  ftheta2 <- ftheta+1
  smoorig2 <- apply(Y, MARGIN = 1, function(x){pForward0(tapes$smotape, ftheta2, x)})
  smopoly2 <- lapply(1:nrow(values$offset), function(i){
    drop(0.5 * ftheta2 %*% matrix(values$Hessian[i, ], ncol = length(ftheta2)) %*% ftheta2 +
      values$offset[i, , drop = FALSE] %*% ftheta2)
  })
  smopoly2 <- unlist(smopoly2)
  expect_equal(smoorig2-smopoly2, constant)
})

