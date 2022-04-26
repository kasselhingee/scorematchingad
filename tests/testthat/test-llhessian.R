test_that("dirichlet ll Hessian matches numerical", {
  beta = c(-0.3, -0.1, 3)
  u <- MCMCpack::rdirichlet(1, beta+1)

  dirichlet_r <- function(u, beta){sum(beta * log(u))}

  lltape <- ptapell(length(beta), length(beta), llname = "dirichlet")
  #Hessian wrt u
  numderiv <- numDeriv::hessian(dirichlet_r, u, beta = beta)
  hess <- pHessian(lltape, u, beta)
  dim(hess) <- c(length(beta), length(beta))
  expect_equal(hess, numderiv, ignore_attr = TRUE)

  #Hessian wrt beta - should be zero
  lltape_theta <- ptapell_theta(lltape, u, beta)
  numderiv <- numDeriv::hessian(function(beta, u){dirichlet_r(u, beta)}, beta, u = u)
  hess <- pHessian(lltape_theta, beta, u)
  dim(hess) <- c(length(beta), length(beta))
  expect_equal(hess, numderiv, ignore_attr = TRUE)
})

test_that("ppi ll Hessian matches numerical", {
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

  ppill_r2 <- function(theta, u){ #theta <- c(diag(ALs), ALs[upper.tri(ALs)], bL, beta0)
    p <- length(u)
    ALs <- matrix(NA, nrow = p - 1, ncol = p - 1)
    diag(ALs) <- theta[1:(p - 1)]
    ALs[upper.tri(ALs)] <- theta[p -1 + 1:((p-2) * (p-1)/2)]
    ALs[lower.tri(ALs)] <- ALs[upper.tri(ALs)]
    bL <- theta[p - 1 + ((p-2) * (p-1)/2) + 1:(p-1)]
    beta0 <- theta[p - 1 + ((p-2) * (p-1)/2) + (p-1) + 1:p]
    ppill_r(u, beta0, ALs, bL)
  }

  stopifnot(isTRUE(all.equal(ppill_r(u, beta0, ALs, bL), ppill_r2(theta, u))))

  lltape <- ptapell(p, length(theta), llname = "ppi")
  #Hessian wrt u
  numderiv <- numDeriv::hessian(ppill_r, u, beta0 = beta0, ALs = ALs, bL = bL)
  hess <- pHessian(lltape, u, theta)
  dim(hess) <- c(p, p)
  expect_equal(hess, numderiv, ignore_attr = TRUE)

  #Hessian wrt beta - should be zero
  lltape_theta <- ptapell_theta(lltape, u, theta)
  numderiv <- numDeriv::hessian(ppill_r2, theta, u = u)
  hess <- pHessian(lltape_theta, theta, u)
  dim(hess) <- c(length(theta), length(theta))
  expect_equal(hess, numderiv, ignore_attr = TRUE)
})
