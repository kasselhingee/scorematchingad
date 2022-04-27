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

ppill_r_S <- function(z, beta0, ALs, bL){
  p=length(z)
  A = matrix(0, nrow = p, ncol = p)
  A[1:p-1, 1:p-1] = ALs
  b = matrix(c(bL, 0), nrow = p, ncol =1)
  out = sum((1 + 2 *beta0) * log(z)) + t(z^2) %*% A %*% (z^2) + t(b) %*% (z^2)
  return(out)
}

# test Jacobian of ll function using numerical differentiation
test_that("ppi likelihood, Jacobian, Hessian for simplex matches numerical estimates", {
  psimplex <- pmanifold("simplex") #because above ppill_r is for the simplex
  lltape <- ptapell(u, theta, llname = "ppi", pman = psimplex)

  # wrt u
  expect_equal(ppill_r(u, beta0, ALs, bL), pForward0(lltape, u, theta), ignore_attr = TRUE)

  #gradiant
  numderiv <- numericDeriv(quote(ppill_r(u, beta0, ALs, bL)), c("u"))
  expect_equal(attr(numderiv,"gradient"), pJacobian(lltape, u, theta), ignore_attr = TRUE, tolerance = 1E-3)

  #Hessian
  numderiv <- numDeriv::hessian(ppill_r, u, beta0 = beta0, ALs = ALs, bL = bL)
  hess <- pHessian(lltape, u, theta)
  dim(hess) <- c(p, p)
  expect_equal(hess, numderiv, ignore_attr = TRUE)

  #wrt u
  lltape_theta <- swapDynamic(lltape, theta, u)
  expect_equal(ppill_r(u, beta0, ALs, bL), pForward0(lltape_theta, theta, u), ignore_attr = TRUE)

  #gradiant
  ppill_r_swap <- function(theta, u){
    pars <- fromPPIparamvec(theta, length(u))
    ppill_r(u, pars$beta0, pars$ALs, pars$bL)
  }
  numderiv <- numericDeriv(quote(ppill_r_swap(theta, u)), c("theta"))
  expect_equal(attr(numderiv,"gradient"), pJacobian(lltape_theta, theta, u),
               ignore_attr = TRUE, tolerance = 1E-7)

  #Hessian
  numderiv <- numDeriv::hessian(ppill_r_swap, theta, u = u)
  hess <- pHessian(lltape_theta, theta, u)
  dim(hess) <- c(length(theta), length(theta))
  expect_equal(hess, numderiv, ignore_attr = TRUE)
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
