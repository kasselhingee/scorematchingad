skip_on_cran() #too slow

##### Prepare First Data Set #####
list2env(ppi_microbiomedata_TCAP(), globalenv())

test_that("hardcoded alr estimator matches cppad calculations for Cyanobacteria/Chloroplast, Actinobacteria, Proteobacteria and pooled", {

# quick check of cppad vs hardcoded
est_hardcoded=ppi(Y = propreal,
         method = "hardcoded", trans = "alr",
         paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0))
hardcodedvals <- ppi_cppad_values(propreal,
         stheta = est_hardcoded$est$paramvec,
         isfixed = t_u2i(ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0)),
         man = "Ralr",
         hsqfun = "ones",
         acut = 1)
expect_lt_v(hardcodedvals$grad, rep(1E-15, length(hardcodedvals$grad)))

skip("next calculation, the cppad estimate, takes hours")
system.time({est_cppad=ppi(Y = propreal,
         method = "cppad", trans = "alr",
         paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
         bdrythreshold = 1E-20, shiftsize = 1E-20,
         control = list(maxit = 1E5, tol = 1E-20 * n))})
# the defaults for Rcgmin are meaning the estimate takes too long to converge!
# from the default starting parameter vec it takes many more iterations than normal to get to converge to the correct result. In this case of default tolerance, then 43961 evaluations of the gradient.
# With current tolerance of 1E-20*n, then 4559 seconds (1.3 hours), 81479 grad evaluations.
expect_equal(est_hardcoded$est$paramvec, est_cppad$est$paramvec, tolerance = 1E-3)
})

test_that("robust ppi via alr estimator matches historical results on dataset with Cyanobacteria/Chloroplast, Actinobacteria, Proteobacteria and pooled", {
# prepare for robustness
#initial values for robust estimators
ALs=matrix(0,p-1,p-1)
bL=matrix(0,p-1,1)
ALs[1,1]= -127480.0929
ALs[1,2]= 14068.39057
ALs[1,3]= 1782.261826
ALs[1,4]=  -240.076568
ALs[2,2]= -8191.17253
ALs[2,3]=  -8.002680
ALs[2,4]= 374.693979
ALs[3,3]= -46.638659
ALs[3,4]= 9.027633
ALs[4,4]= -39.208915
ALs[lower.tri(ALs)] <- t(ALs)[lower.tri(ALs)]
beta0=matrix(0,p,1)
beta0[1]=-0.80
beta0[2]=-0.85
beta0[3]=0
beta0[4]=-0.2
beta0[5]=0
bL_est=bL
ALs_est=ALs
beta0_est=beta0
sp=p-1

#try simulating to see required maxden
expect_silent(sim <- rppi(1000, beta0, ALs, bL, maxden = 0))


#calculate robust estimates
cW=0.7
cWvec <- ppi_cW(cW, TRUE, TRUE, FALSE, FALSE, FALSE)
est1=ppi_robust(Y = propreal,
                cW = cWvec,
                method = "hardcoded", trans = "alr",
                paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                paramvec_start = ppi_paramvec(AL = ALs_est, bL = bL_est, beta = beta0_est))
#estimate of A_L:
expect_snapshot_value(signif(est1$est$ALs,6), style = "json2")
#estimate of beta:
expect_snapshot_value(signif(est1$est$beta,6), style = "json2")


# check that restarting fp ends up in a quick finish
est1b=ppi_robust(Y = propreal,
                cW = cWvec,
                method = "hardcoded", trans = "alr",
                paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                paramvec_start = est1$est$paramvec,
                fpcontrol = list(ConvergenceMetricThreshold = 1E-9)) #slightly larger to beat est1
expect_equal(est1b$est$paramvec, est1$est$paramvec)

# Doing the full fp search with cppad fitting took too long and with maxit = 100, the estimate was poor.
# instead verify the est1 result. The Windham solution theta is such that tauc(theta) = T(F weighted by c,theta)
ldenfun <- function(Y, theta){ #here theta is the usual parameters of PPI model from
  return(drop(dppi(Y, paramvec = theta)))
}

weights <- Windham_weights(ldenfun = ldenfun, Y = propreal,
                theta = est1$est$paramvec,
                cW = cWvec)
expect_lt(max(abs(est1$info$finalweights - weights)), 1E-10)

vals <- ppi_cppad_values(prop = propreal,
         stheta = est1$est$paramvec * (1 + cWvec), #this is Tauc(theta), which should be the solution to the weighted estimator
         isfixed = t_u2i(ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0)),
         man = "Ralr",
         hsqfun = "ones",
         acut = 1,
         w = weights)
expect_lt(sum(vals$grad^2), 1E-10)
})

#### Test Second Data Set ####

test_that("robust ppi via alr estimator matches historical results on dataset with Spirochates, Verrucomicrobia, Cyanobacteria/Chloroplast, TM7 and pooled", {
  data("microdata", package = "scorecompdir")
  countdata=as.matrix(microdata[,12:31])

  #sample size
  n=94

  #dimension
  p=20

  #calculate totals
  tot=matrix(0,n,1)
  for (j in 1:p)
  {
    tot=tot+countdata[,j]
  }
  tot=as.vector(tot)

  #proportion data
  prop=countdata
  for (j in 1:n)
  {
    prop[j,]=countdata[j,]/tot[j]
  }

  ##Reduce dimensions to p=5


  #dimension
  p=5

  #calculate 5D dataset
  comb=matrix(0,n,p)
  comb[,1]=prop[,14]
  comb[,2]=prop[,18]
  comb[,3]=prop[,5]
  comb[,4]=prop[,16]
  for (j in 1:sum(p,-1))
  {
    comb[,p]=comb[,p]+comb[,j]
  }
  comb[,p]=1-comb[,p]

  #save data
  propreal=comb
  proprealA=propreal

  ##Estimation

  #initial values for robust estimators
  ALs=matrix(-10000,p-1,p-1)
  bL=matrix(0,p-1,1)
  beta0=matrix(0,p,1)
  beta0[1]=-0.80
  beta0[2]=-0.80
  beta0[3]=-0.80
  beta0[4]=-0.80
  bL_est=bL
  ALs_est=ALs
  beta0_est=beta0

  #calculate robust estimates
  cW=1.25
  est1=ppi_robust(Y = propreal,
                   cW = ppi_cW(cW, TRUE, TRUE, TRUE, TRUE, FALSE),
                   method = "hardcoded", trans = "alr",
                   paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                   paramvec_start = ppi_paramvec(AL = ALs_est, bL = bL_est, beta = beta0_est))
  #estimate of A_L:
  expect_snapshot_value(signif(fromPPIparamvec(est1$est$paramvec)$ALs,6), style = "json2")
  #estimate of beta:
  expect_snapshot_value(signif(fromPPIparamvec(est1$est$paramvec)$beta,6), style = "json2")
})

