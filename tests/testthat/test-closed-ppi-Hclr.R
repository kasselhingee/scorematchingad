
test_that("ppi() warns when bdrythreshold too low", {
  set.seed(12345)
  model <- ppi_egmodel(100, maxden = 4)

  expect_warning(ppi(Y = model$sample, trans = "clr", bdrythreshold = 1E-10),
                 "bdrythreshold")
  expect_warning(ppi_smvalues(Y = model$sample, evalparam = model$theta, trans = "clr", bdrythreshold = 1E-10),
                 "bdrythreshold")
})

test_that("Fitting ppi via clr transform gets close estimates via sqrt", {
  set.seed(1234)
  model <- ppi_egmodel(100, maxden = 4)

  out <- ppi(Y = model$sample,
             paramvec = ppi_paramvec(p = 3, betap = -0.5),
             trans = "clr", bdrythreshold = 1E-5)
  expect_absdiff_lte_v(out$est$paramvec, model$theta, out$SE$paramvec * 3)
  out2 <- ppi(Y = model$sample,
             paramvec = ppi_paramvec(p = 3, betap = -0.5),
             trans = "sqrt", acut = 0.01, divweight = "minsq")
  expect_absdiff_lte_v(out$est$paramvec, out2$est$paramvec, 
       out$SE$paramvec * c(1, 0.1, 0.5, 1, 0.1, 1, 1, 1))
})


test_that("Fitting ppi all parameters via clr transform can get close to true values for positive-looking beta", {
  skip_on_cran()
  set.seed(12345)
  model <- ppi_egmodel(10000, maxden = 4)

  out <- ppi(Y = model$sample,
             trans = "clr", bdrythreshold = 1E-5)

  expect_absdiff_lte_v(out$est$paramvec, model$theta, out$SE$paramvec * 3)
  # and that the SE are small
  expect_lt(max(abs(out$SE$paramvec), na.rm = TRUE), 20)
})
  # method of moments from JASA paper supplementary
  dirichmom <- function(X) {
    # Method of Moments estimates of the Dirichlet Distribution
    temp <- dim(X); n <- temp[1]; m <- temp[2]
    # X <- cbind(X,matrix(1-apply(X,1,sum)))
    mom <- apply(X,2,mean)*(mean(X[,1])-mean(X[,1]^2))/(mean(X[,1]^2) - ((mean(X[,1]))^2))
    return(mom - 1)
  }

test_that("Fitting Dirichlet with no bdrythreshold fails for high concentrations", {
  set.seed(13415)
  beta0 <- c(-0.7, -0.1, 0) #the largest element is going to look a bit like beta>0 
  Y <- MCMCpack::rdirichlet(10000, beta0+1)
  suppressWarnings({out <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0, betap = tail(beta0, 1)),           
             trans = "clr", bdrythreshold = 0)})

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
  expect_error(suppressWarnings({out <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0, betap = tail(beta0, 1)),
             trans = "clr", bdrythreshold = 0)}))
 
  out <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0, betap = tail(beta0, 1)),
             trans = "clr", bdrythreshold = 0.1)
})


test_that("Fitting Dirichlet with proper bdrythreshold gets close to true values and other estimators for high concentrations", {
  set.seed(13412)
  beta0 <- c(-0.5, 0, -0.8) #the largest element is going to look a bit like beta>0
  Y <- MCMCpack::rdirichlet(10000, beta0+1)
  out <- ppi(Y = Y,
             paramvec = ppi_paramvec(p = 3, AL=0, bL = 0),
             trans = "clr", bdrythreshold = 1E-5)
  out$est$beta
 
  expect_absdiff_lte_v(out$est$beta, beta0, out$SE$beta * 3)
  expect_lt(max(abs(out$SE$beta), na.rm = TRUE), 1E-1)

  # and result close to dirichlet moment estimator
  expect_absdiff_lte_v(out$est$beta, dirichmom(Y), c(0.01, 0.02, 0.01))
})


test_that("Hess + Offset match gradient for Hclr in interior", {
  set.seed(1245)
  mod <- ppi_egmodel(100)
  Y <- mod$sample
  Y <- rbind(Y, rep(1/3, 3))

  tapes <- buildsmotape("Hclr", "ppi", utape = c(0.2, 0.3, 0.5), 
                        usertheta = ppi_paramvec(p = 3))

  # find boundary points and remove them - expect results pJacobian to give good results when no points need approximation
  isbdry <- simplex_isboundary(Y, 1E-5)
  Y <- Y[!isbdry, ]
  values <- quadratictape_parts(tapes$smotape, Y)

  # expect results to match for gradient
  gradorig <- t(apply(Y, MARGIN = 1, function(x){pJacobian(tapes$smotape, mod$theta, x)}))

  gradpoly <- lapply(1:nrow(values$offset), function(i){
    drop(matrix(values$Hessian[i, ], ncol = length(mod$theta)) %*% mod$theta + 
      t(values$offset[i, , drop = FALSE]))
  })
  gradpoly <- do.call(rbind, gradpoly)
  expect_equal(gradorig, gradpoly, tolerance = 1E-5)
})

test_that("Taylor approximations accuracy study", {
  set.seed(1245)
  mod <- ppi_egmodel(1000)
  Y <- mod$sample
  tapes <- buildsmotape("Hclr", "ppi", utape = c(0.2, 0.3, 0.5), 
                        usertheta = ppi_paramvec(p = 3))

  # find boundary points 
  isbdry <- simplex_isboundary(Y, 1E-3)
  Yapproxcentres <- Y
  Yapproxcentres[!isbdry, ] <- NA
  Yapproxcentres[isbdry, ] <- simplex_boundaryshift(Y[isbdry, , drop = FALSE], shiftsize = 1E-5)

  approxedvals <- tape_smvalues(tapes$smotape, Y, mod$theta, xcentres = Yapproxcentres, approxorder = 10)
  naivevals <- tape_smvalues(tapes$smotape, Y, mod$theta)

  #compute distance to boundary
  bdrydist <- apply(Y, MARGIN = 1, min)
  bdrycomponent <- apply(Y, MARGIN = 1, which.min)

  # look at obj values
  library(dplyr)
  library(ggplot2)
  data.frame(naive = naivevals$obj,
             approxed = approxedvals$obj,
             bdrydist = bdrydist) %>%
    ggplot(aes(x = bdrydist, y = abs(naive - approxed))) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() + 
    geom_smooth(method = "lm", se = FALSE)
  # result is straight in log-log (i.e. power law)
  fit <- lm(log(err) ~ log(bdrydist), data = data.frame(err = abs(naivevals$obj - approxedvals$obj)[isbdry, ], bdrydist = bdrydist[isbdry])) 
  # about is err ~ exp(-40) / bdrydist^2
  # reaches 1E-6 in error just below 1E-5 in bdry distance

  # look at grad values
  data.frame(err = abs(naivevals$grad-approxedvals$grad),
           bdrydist = bdrydist,
           bdrycomponent = factor(bdrycomponent)) %>%
    tidyr::pivot_longer(!c(bdrydist, bdrycomponent)) %>%
    ggplot(aes(x = bdrydist, y = value, col = bdrycomponent)) +
    facet_wrap(vars(name)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = "lm", se = FALSE)
  # again straight lines in the log-log plots
  # most components of the grad reach error of 1E-13 at about bdry dist of 1E-05  # the beta components are at this limit at 1E-3 bdry distance.
  # I suspect the naive calculations are the one that are wrong here (since error is consistent relation to bdry distance) and that the cause is some numerical instability

  # look at hess values
  data.frame(err = abs(naivevals$hess-approxedvals$hess),
           bdrydist = bdrydist,
           bdrycomponent = factor(bdrycomponent)) %>%
    tidyr::pivot_longer(!c(bdrydist, bdrycomponent)) %>%
    ggplot(aes(x = bdrydist, y = value, col = bdrycomponent)) +
    facet_wrap(vars(name)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10()
  # relationship between accuracy and distance (linear in log-log) is a power law for the errors above 1E-14.
})

test_that("W is symmetric for ppi with clr, fitting all parameters", {
  # if W is symmetric then should equal the smo up to a constant wrt theta
  # and my expectation is that W is symmetric for ppi with Hclr
  set.seed(1245)
  mod <- ppi_egmodel(100)
  Y <- mod$sample
  Y <- rbind(Y, rep(1/3, 3))
  isbdry <- simplex_isboundary(Y, 1E-5)
  Y <- Y[!isbdry, ]

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


test_that("printgraph() runs", {
  set.seed(1245)
  beta0 <- c(-0.7, -0.1, -0.5, -0.1, -0.6) 
  Y <- MCMCpack::rdirichlet(10, beta0+1)
  tapes <- buildsmotape("Hclr", "ppi", utape = Y[1, ], 
                        usertheta = ppi_paramvec(p = 5))

  expect_snapshot(printgraph(tapes$smotape))
  Jtape <- pTapeJacobian(tapes$smotape, ppi_paramvec(AL = 0, bL = 0, beta = beta0), Y[1, ])
  expect_snapshot(printgraph(Jtape))
  Htape <- pTapeHessian(tapes$smotape, ppi_paramvec(AL = 0, bL = 0, beta = beta0), Y[1, ])
  expect_snapshot(printgraph(Htape))
})
