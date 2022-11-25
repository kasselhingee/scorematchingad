skip_on_cran() #too slow

##### Prepare First Data Set #####
list2env(ppi_microbiomedata_TCAP(), globalenv())

test_that("hardcoded alr estimator matches cppad calculations for Cyanobacteria/Chloroplast, Actinobacteria, Proteobacteria and pooled", {

# quick check of closed vs hardcoded
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

  est_cppad=ppi(Y = propreal,
         method = "closed", trans = "alr",
         paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
         bdrythreshold = 1E-20, shiftsize = 1E-20,
         control = list(maxit = 1E5, tol = 1E-20 * n))
  expect_equal(est_hardcoded$est$paramvec, est_cppad$est$paramvec, ignore_attr = TRUE)
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

# check that ppi_robust with 'closed' method also works
est1c=ppi_robust(Y = propreal,
                cW = cWvec,
                method = "closed", trans = "alr",
                paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                paramvec_start =  ppi_paramvec(AL = ALs_est, bL = bL_est, beta = beta0_est)) 
  expect_equal(est1c$est$paramvec, est1$est$paramvec)
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

  # check that closed method works too
  est1b=ppi_robust(Y = propreal,
                   cW = ppi_cW(cW, TRUE, TRUE, TRUE, TRUE, FALSE),
                   method = "closed", trans = "alr",
                   paramvec = ppi_paramvec(p=ncol(propreal), bL = 0, betap = 0),
                   paramvec_start = ppi_paramvec(AL = ALs_est, bL = bL_est, beta = beta0_est))
  expect_equal(est1b$est$paramvec, est1$est$paramvec)
})

