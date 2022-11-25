#### Preparing microbiome data #############################################

list2env(ppi_microbiomedata_cleaned_TCAP(), globalenv())

#a_c for h function:
acut=0.01

#### Test including b_L ####
test_that("estimator1 and SE is historically correct with b_L included (article Table 4)", {


  #calculate scoring estimate:
  estimator= estimator1(propreal,acut,1, beta0, computeSE = TRUE)
  estimate1=estimator$est$paramvec
  #rearrange to historical ordering
  ordindx <- order(combparam2uppertriorder(length(estimate1))) #the plus p is for the beta that isn't estimated
  estimate1 <- estimate1[ordindx][1:(length(estimate1) - p)]
  dim(estimate1) <- c(length(estimate1), 1)
  #check historically
  expect_snapshot_value(signif(estimate1, 8), style = "json2") #8 is the default number of digits for jsonlite::serializeJSON

  #check it matches cppad ppi()
  est_cppad <- ppi(Y = propreal, acut = acut,
                   method = "closed",
                   trans = "sqrt", divweight = "minsq",
                   bdrythreshold = 1E-15, shiftsize = 1E-10,
                   paramvec = ppi_paramvec(beta = beta0),
                   control = list(maxit = 100000, tol = 1E-20))
  expect_equal(est_cppad$est$paramvec, estimator$est$paramvec, tolerance = 1E-3)

  #estimate of W matrix
  W_est=estimator$info$W
  expect_snapshot_value(round(max(W_est), 8), style = "json2") #have to use round here because the json conversion doesn't necessarily show it in scientific notation
  expect_snapshot_value(signif(mean(W_est), 8), style = "json2")
  expect_snapshot_value(signif(which.max(W_est), 8), style = "json")

  #standard errors
  std1= estimator$SE$paramvec
  #rearrange back to combn ordering
  std1 <- std1[ordindx][1:(length(std1) - p)]
  expect_snapshot_value(signif(std1, 8), style = "json2")

  #estimated parameters
  thetamats <- fromPPIparamvec(estimator$est$paramvec)
  ALs <- thetamats$ALs
  bL <- thetamats$bL
  dim(bL) <- c(length(bL), 1)
  #values in Table 4 in the article:
  expect_snapshot_value(signif(ALs[upper.tri(ALs, diag = TRUE)], 8), style = "json2")
  expect_snapshot_value(signif(bL, 8), style = "json2")
  expect_snapshot_value(signif(estimate1/std1, 8), style = "json2")
})


test_that("alr and cppad estimator for this data set are consistent", {
  #check alr estimators too
  est_alr <- ppi(Y = propreal, method = "hardcoded",
                 trans = "alr", 
                 paramvec = ppi_paramvec(p = ncol(propreal), bL = 0, betap = tail(beta0, 1)))

  skip("next calculation, the cppad estimate, takes hours")
  stop("The number iterations required has not been checked - has taken multiple hours without finishing")
  system.time({est_alr_cppad <- ppi(Y = propreal, method = "closed",
                 trans = "alr", 
                 bdrythreshold = 1E-15, shiftsize = 1E-15,
                 approxorder = 10,
                 control = list(maxit = 1E15, tol = 1E-20 * nrow(propreal)),
                 paramvec = ppi_paramvec(p = ncol(propreal), bL = 0, betap = tail(beta0, 1)))})
  expect_equal(est_alr$est$paramvec, est_alr_cppad$est$paramvec, tolerance = 1E-2)


})

#### Test omitting b_L ####
test_that("estimator1 and SE is historically correct with b_L ommitted (article table 3)", {

  #calculate scoring estimate:
  estimator=estimator1(propreal,acut,0, beta0, computeSE = TRUE)
  estimate1=estimator$est$paramvec
  #rearrange to historical ordering
  ordindx <- order(combparam2uppertriorder(length(estimate1)))
  estimate1 <- estimate1[ordindx][1:(length(estimate1) - p - (p-1))]
  dim(estimate1) <- c(length(estimate1), 1)
  expect_snapshot_value(signif(estimate1, 8), style = "json2")

  #estimate of W matrix
  W_est=estimator$info$W
  expect_snapshot_value(round(max(W_est), 8), style = "json2") #have to use round here because the json conversion doesn't necessarily show it in scientific notation
  expect_snapshot_value(signif(mean(W_est), 8), style = "json2")
  expect_snapshot_value(signif(which.max(W_est), 8), style = "json")

  #standard errors
  std1= estimator$SE$paramvec
  std1 <- std1[ordindx][1:(length(estimate1))] #rearrange to combn ordering for historical comparison

  #estimated parameters
  thetamats <- fromPPIparamvec(estimator$est$paramvec)
  ALs <- thetamats$ALs

  #values in Table 3 in the article:
  expect_snapshot_value(signif(ALs[upper.tri(ALs, diag = TRUE)], 8), style = "json2")
  expect_snapshot_value(round(estimate1/std1, 8), style = "json2")

})

