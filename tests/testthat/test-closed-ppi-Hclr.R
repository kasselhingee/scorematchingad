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


test_that("Fitting Dirichlet via clr transform gets close to true values and other estimators", {
  skip_on_cran()
  set.seed(13412)
  beta0 <- c(-0.5, 0, -0.8) #the largest element is going to look a bit like beta>0
  Y <- MCMCpack::rdirichlet(10000, beta0+1)
  out <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0),
             trans = "clr", bdrythreshold = 0)
  out$est$beta
 
  library(ggtern)
  library(ggplot2)
  library(dplyr)
  colnames(Y) <- c("x", "y", "z")
  Y %>%
   as.data.frame() %>%
   ggplot(aes(x=x, y = y, z = z)) +
   coord_tern() +
   geom_density_tern(bdl = 0.1, bdl.val = NA) +
   geom_point(shape = "+")

  clr(Y) %>%
   as.data.frame() %>%
   ggplot() +
   geom_point(aes(x = x, y = y))

  
  expect_absdiff_lte_v(out$est$beta, beta0, out$SE$beta * 3)
  # but the estimate is wrong by 10 orders of magnitude
  expect_lt(max(abs(out$SE$beta), na.rm = TRUE), 20)

  # vs method of moments from JASA paper supplementary
  dirichmom <- function(X) {
    # Method of Moments estimates of the Dirichlet Distribution
    temp <- dim(X); n <- temp[1]; m <- temp[2]
    # X <- cbind(X,matrix(1-apply(X,1,sum)))
    mom <- apply(X,2,mean)*(mean(X[,1])-mean(X[,1]^2))/(mean(X[,1]^2) - ((mean(X[,1]))^2))
    return(matrix(mom))
  }
  dirichmom(Y)-1 #this is very accurate! :)

  out2 <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0),
             trans = "sqrt", divweight = "minsq", acut = 0.01)

  out2$est$beta

})

test_that("Fitting ppi via clr transform with fixed beta gets close to true values and other estimators", {
  skip_on_cran()
  set.seed(1234)
  model <- ppi_egmodel(100000, maxden = 4)

  out <- ppi(Y = model$sample,
             paramvec = ppi_paramvec(beta = model$beta),
             trans = "clr")
  expect_absdiff_lte_v(out$est$paramvec, model$theta, out$SE$paramvec * 3)
  # but the estimate is wrong by 10 orders of magnitude
  expect_lt(max(abs(out$SE$paramvec), na.rm = TRUE), 20)


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

test_that("tapefromM and toM matches R-computed versions", {
  set.seed(1245)
  mod <- ppi_egmodel(100)
  Y <- mod$sample
  Y <- rbind(Y, rep(1/3, 3))

  Hclr <- manifoldtransform("Hclr")

  # check evaluation
  clrY <- clr(Y)
  expect_equal(clrY, t(apply(Y, MARGIN = 1, function(u){ptoM(Hclr, u)})))

  tapefromM <- ptapefromM(clrY[1, ], Hclr)
  vals <- tape_eval(tapefromM, clrY, matrix(nrow = nrow(Y), ncol = 0))
  expect_equal(vals, Y)
})

test_that("tapefromM Jacobian and taped logdetJfromM match R-computed versions", {
  set.seed(1245)
  mod <- ppi_egmodel(100)
  Y <- mod$sample
  Y <- rbind(Y, rep(1/3, 3))

  Hclr <- manifoldtransform("Hclr")

  # Start of check Jacobian
  fromMJ <- t(apply(clrY, MARGIN = 1, function(x){pJacobian(tapefromM, x, matrix(nrow = 0, ncol = 0))}))
  fromMJnum <- t(apply(clrY, MARGIN = 1, function(x){
    b <- x
    J <- numericDeriv(quote(pForward0(tapefromM, b , matrix(nrow = 0, ncol = 0))), "b")
    return(as.vector(attr(J, "gradient")))
}))
  expect_equal(fromMJ, fromMJnum, tolerance = 1E-6)


  # Check the determinants: via det of Jacobian above, via numerical via R analytical, and via CppAD tape
  # determinant calculated from taped Jacobian
  detJ <- apply(fromMJ, MARGIN = 1, function(J){
    Jmat <- matrix(J, nrow = sqrt(length(J)), ncol = sqrt(length(J)))
    return(det(Jmat))
  })
  detJfromMtape <- ptapelogdetJ(tapefromM, clrY[1,], vector(mode = "numeric", length = 0))
  # log determinant evaluation via tape
  ldetJ_cppad <- tape_eval(detJfromMtape, clrY, matrix(nrow = nrow(Y), ncol = 0))
  expect_equal(ldetJ_cppad, log(abs(detJ)), ignore_attr = TRUE, tolerance = 1E-7)

  # analytic determinant
  ldetJ_Rdirect <- apply(clrY, MARGIN = 1, ldetJfromM)
  expect_equal(log(abs(detJ)), ldetJ_Rdirect, tolerance = 1E-6)
  expect_equal(ldetJ_cppad, ldetJ_Rdirect, ignore_attr = TRUE, tolerance = 1E-6)
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

