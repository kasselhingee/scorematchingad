clr <- function(Y){
  logY <- log(Y)
  lgeommean <- rowSums(logY)/3
  return(logY - lgeommean)
}

  # vs method of moments from JASA paper supplementary
  dirichmom <- function(X) {
    # Method of Moments estimates of the Dirichlet Distribution
    temp <- dim(X); n <- temp[1]; m <- temp[2]
    # X <- cbind(X,matrix(1-apply(X,1,sum)))
    mom <- apply(X,2,mean)*(mean(X[,1])-mean(X[,1]^2))/(mean(X[,1]^2) - ((mean(X[,1]))^2))
    return(mom - 1)
  }

test_that("Fitting Dirichlet via clr transform with fixed betap gets close to true values and other estimators for low concentrations", {
  skip_on_cran()
  set.seed(13412)
  beta0 <- c(-0.5, -0.1, 0) #the largest element is going to look a bit like beta>0
  Y <- MCMCpack::rdirichlet(10000, beta0+1)
  out <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0, betap = tail(beta0, 1)),
             trans = "clr", bdrythreshold = 0)
  
  expect_absdiff_lte_v(out$est$beta, beta0, out$SE$beta * 3)
  expect_lt(max(abs(out$SE$beta), na.rm = TRUE), 1)

  expect_equal(dirichmom(Y), out$est$beta, tolerance = 0.1)
  out2 <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0),
             trans = "sqrt", divweight = "minsq", acut = 0.01)

  out2$est$beta
})

test_that("Fitting Dirichlet via clr transform with fixed betap fails for high concentrations", {
  skip_on_cran()
  set.seed(13415)
  beta0 <- c(-0.7, -0.1, 0) #the largest element is going to look a bit like beta>0 
  Y <- MCMCpack::rdirichlet(10000, beta0+1)
  out <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0, betap = tail(beta0, 1)),           
             trans = "clr", bdrythreshold = 0)

  # estimate is very poor
  expect_gt(max(abs(dirichmom(Y) - out$est$beta)), 0.9)
  # standard errors fail a lot
  expect_error(expect_absdiff_lte_v(out$est$beta, beta0, out$SE$beta * 10))
  out2 <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0),
             trans = "sqrt", divweight = "minsq", acut = 0.01)
  #sqrt estimator performs well still, even with unfixed beta
  expect_absdiff_lte_v(out2$est$beta, beta0, out2$SE$beta * 3)
  expect_lt(max(abs(out2$SE$beta), na.rm = TRUE), 0.2)

  # even higher concentration has an invertibility problem
  set.seed(1341)
  beta0 <- c(-0.8, -0.1, 0) #the largest element is going to look a bit like beta>0 
  Y <- MCMCpack::rdirichlet(10000, beta0+1)
  expect_error({out <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0, betap = tail(beta0, 1)),
             trans = "clr", bdrythreshold = 0)})
 
  out <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0, betap = tail(beta0, 1)),
             trans = "clr", bdrythreshold = 0.1)
})


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

test_that("Fitting ppi via clr transform with fixed beta gets close to true values", {
  skip_on_cran()
  set.seed(1234)
  model <- ppi_egmodel(100000, maxden = 4)

  out <- ppi(Y = model$sample,
             paramvec = ppi_paramvec(beta = model$beta),
             trans = "clr")
  expect_absdiff_lte_v(out$est$paramvec, model$theta, out$SE$paramvec * 3)
  # and that the SE are small
  expect_lt(max(abs(out$SE$paramvec), na.rm = TRUE), 20)
})


test_that("Fitting ppi via clr transform with unfixed beta gets close to true values", {
  skip_on_cran()
  set.seed(12345)
  model <- ppi_egmodel(1000000, maxden = 4)

  out <- ppi(Y = model$sample,
             trans = "clr")

  expect_absdiff_lte_v(out$est$paramvec, model$theta, out$SE$paramvec * 3)
  # and that the SE are small
  expect_lt(max(abs(out$SE$paramvec), na.rm = TRUE), 20)
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

