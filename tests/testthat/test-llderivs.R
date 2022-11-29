set.seed(1234) #chosen so that u is far from boundary - needed for numericDeriv::Hessian
m = ppi_egmodel_p4(1, 8)
p = m$p
u = m$sample
ALs = m$ALs
bL = m$bL
beta0 = m$beta0
theta = m$theta

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
  out = sum((1 + 2 *beta0) * log(z)) + t(z^2) %*% A %*% (z^2) + t(b) %*% (z^2) + log(2) * p
  return(out)
}

# test Jacobian of ll function using numerical differentiation
test_that("ppi likelihood, Jacobian, Hessian for simplex matches numerical estimates", {
  psimplex <- manifoldtransform("simplex") #because above ppill_r is for the simplex
  lltape <- ptapell(u, theta, llname = "ppi", pman = psimplex, fixedtheta = rep(FALSE, length(theta)), verbose = FALSE)

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
    ppill_r(u, pars$beta, pars$ALs, pars$bL)
  }
  numderiv <- numericDeriv(quote(ppill_r_swap(theta, u)), c("theta"))
  expect_equal(attr(numderiv,"gradient"), pJacobian(lltape_theta, theta, u),
               ignore_attr = TRUE, tolerance = 1E-7)

  #Hessian
  numderiv <- numDeriv::hessian(ppill_r_swap, theta, u = u)
  hess <- pHessian(lltape_theta, theta, u)
  dim(hess) <- c(length(theta), length(theta))
  expect_equal(hess, numderiv, ignore_attr = TRUE, tolerance = 1E-5)
})

# test Jacobian of ll function using numerical differentiation
test_that("ppi likelihood, Jacobian, Hessian for sphere matches numerical estimates", {
  psphere <- manifoldtransform("sphere")
  lltape <- ptapell(u, theta, llname = "ppi", pman = psphere, fixedtheta = rep(FALSE, length(theta)), verbose = FALSE)

  # wrt u
  expect_equal(ppill_r_S(u, beta0, ALs, bL), pForward0(lltape, u, theta), ignore_attr = TRUE)

  #gradiant
  numderiv <- numericDeriv(quote(ppill_r_S(u, beta0, ALs, bL)), c("u"))
  expect_equal(attr(numderiv,"gradient"), pJacobian(lltape, u, theta), ignore_attr = TRUE, tolerance = 1E-3)

  #Hessian
  numderiv <- numDeriv::hessian(ppill_r_S, u, beta0 = beta0, ALs = ALs, bL = bL)
  hess <- pHessian(lltape, u, theta)
  dim(hess) <- c(p, p)
  expect_equal(hess, numderiv, ignore_attr = TRUE)

  #wrt theta
  lltape_theta <- swapDynamic(lltape, theta, u)
  expect_equal(ppill_r_S(u, beta0, ALs, bL), pForward0(lltape_theta, theta, u), ignore_attr = TRUE)

  #gradiant
  ppill_r_S_swap <- function(theta, z){
    pars <- fromPPIparamvec(theta, length(z))
    ppill_r_S(z, pars$beta, pars$ALs, pars$bL)
  }
  numderiv <- numericDeriv(quote(ppill_r_S_swap(theta, u)), c("theta"))
  expect_equal(attr(numderiv,"gradient"), pJacobian(lltape_theta, theta, u),
               ignore_attr = TRUE, tolerance = 1E-7)

  #Hessian
  numderiv <- numDeriv::hessian(ppill_r_S_swap, theta, z = u)
  hess <- pHessian(lltape_theta, theta, u)
  dim(hess) <- c(length(theta), length(theta))
  expect_equal(hess, numderiv, ignore_attr = TRUE, tolerance = 1E-5)
})


test_that("dirichlet ll evaluation and Jacobian matches expected", {
  beta = c(-0.3, -0.1, 3)
  set.seed(1234)
  u <- MCMCpack::rdirichlet(1, beta+1)

  dirichlet_r <- function(u, beta){sum(beta * log(u))}

  psimplex <- manifoldtransform("simplex")
  lltape <- ptapell(u, beta, llname = "dirichlet", pman = psimplex, fixedtheta = rep(FALSE, length(beta)), verbose = FALSE)
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
