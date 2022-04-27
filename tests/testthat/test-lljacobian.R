# test Jacobian of ll function using numerical differentiation
test_that("ppi likelihood Jacobian matches numerical estimates", {
  p=4
  #parameters for the PPI model
  muL=matrix(0,p-1,1)
  muL[1:sum(p,-1)]=0.12
  aa=matrix(1/500,p-1,1)
  D=diag(as.vector(aa))
  SigA=D
  SigA[1,1]=SigA[1,1]*2
  cor=0.5
  SigA[1,2]=cor*sqrt(SigA[1,1]*SigA[2,2])
  SigA[2,1]=SigA[1,2]
  ALs=-0.5*solve(SigA)
  bL=solve(SigA)%*%muL
  beta0=matrix(-0.8,p,1)
  beta0[p]=-0.5
  theta <- c(diag(ALs), ALs[upper.tri(ALs)], bL, beta0)

  set.seed(123)
  u <- cdabyppi:::rhybrid(1,p,beta0,ALs,bL,4)$samp3

  ppill_r <- function(u, beta0, ALs, bL){
    p=length(u)
    A = matrix(0, nrow = p, ncol = p)
    A[1:p-1, 1:p-1] = ALs
    b = matrix(c(bL, 0), nrow = p, ncol =1)
    out = sum(beta0 * log(u)) + t(u) %*% A %*% u + t(b) %*% u
    return(out)
  }

  psimplex <- pmanifold("simplex") #because above ppill_r is for the simplex
  lltape <- ptapell(u, theta, llname = "ppi", pman = psimplex)
  expect_equal(ppill_r(u, beta0, ALs, bL), pForward0(lltape, u, theta), ignore_attr = TRUE)

  #gradiant
  numderiv <- numericDeriv(quote(ppill_r(u, beta0, ALs, bL)), c("u"))
  expect_equal(attr(numderiv,"gradient"), pJacobian(lltape, u, theta), ignore_attr = TRUE, tolerance = 1E-3)

  #gradient wrt theta
  lltape_theta <- swapDynamic(lltape, theta, u)
  numderiv <- numericDeriv(quote(ppill_r(u,  beta0, ALs, bL)), c("ALs", "bL", "beta0"))
  numgrad = attr(numderiv,"gradient")
  inthetaform = c(numgrad[c(1, 5, 9)], #diag of ALs
                  sum(numgrad[c(2, 4)]), #1,2 and 2,1 element of ALs
                  sum(numgrad[c(3, 7)]), #1,3 and 3,1 element of ALs
                  sum(numgrad[c(6, 8)]), #2,3 and 3,2 element of ALs
                  numgrad[10:12], #bL
                  numgrad[13:16]) #beta0
  expect_equal(inthetaform, pJacobian(lltape_theta, theta, u),
               ignore_attr = TRUE, tolerance = 1E-7)
})

test_that("dirichlet ll evaluation and Jacobian matches expected", {
  beta = c(-0.3, -0.1, 3)
  u <- MCMCpack::rdirichlet(1, beta+1)

  dirichlet_r <- function(u, beta){sum(beta * log(u))}

  psimplex <- pmanifold("simplex")
  lltape <- ptapell(u, beta, llname = "dirichlet", pman = psimplex)
  #forward0
  expect_equal(dirichlet_r(u, beta), pForward0(lltape, u, beta), ignore_attr = TRUE)

  #gradiant wrt u
  numderiv <- numericDeriv(quote(dirichlet_r(u, beta)), c("u"))
  expect_equal(attr(numderiv,"gradient"), pJacobian(lltape, u, beta), ignore_attr = TRUE)

  #gradient wrt beta
  lltape_theta <- swapDynamic(lltape, beta, u)
  numderiv <- numericDeriv(quote(dirichlet_r(u, beta)), c("beta"))
  expect_equal(attr(numderiv,"gradient"), pJacobian(lltape_theta, beta, u), ignore_attr = TRUE)
})
