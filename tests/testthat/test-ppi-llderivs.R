set.seed(1234) #chosen so that u is far from boundary - needed for numericDeriv::Hessian
m = rppi_egmodel_p4(1, 8)
p = m$p
u = as.vector(m$sample)
AL = m$AL
bL = m$bL
beta0 = m$beta
theta = m$theta

ppill_r <- function(u, beta0, AL, bL){
  p=length(u)
  A = matrix(0, nrow = p, ncol = p)
  A[1:p-1, 1:p-1] = AL
  b = matrix(c(bL, 0), nrow = p, ncol =1)
  out = sum(beta0 * log(u)) + t(u) %*% A %*% u + t(b) %*% u
  return(out)
}

ppill_r_S <- function(z, beta0, AL, bL){
  p=length(z)
  A = matrix(0, nrow = p, ncol = p)
  A[1:p-1, 1:p-1] = AL
  b = matrix(c(bL, 0), nrow = p, ncol =1)
  out = sum((1 + 2 *beta0) * log(z)) + t(z^2) %*% A %*% (z^2) + t(b) %*% (z^2) + log(2) * p
  return(out)
}

# test Jacobian of ll function using numerical differentiation
test_that("ppi likelihood, Jacobian, Hessian for simplex matches numerical estimates", {
  psimplex <- manifoldtransform("sim", "identity", "sim") #because above ppill_r is for the simplex
  lltape <- tapell("ppi", u, NA * theta, tran = psimplex$tran, function(n){theta})

  # wrt u
  expect_equal(ppill_r(u, beta0, AL, bL), lltape$eval(u, theta), ignore_attr = TRUE)

  #gradiant
  numderiv <- numericDeriv(quote(ppill_r(u, beta0, AL, bL)), c("u"))
  expect_equal(attr(numderiv,"gradient"), lltape$Jac(u, theta), ignore_attr = TRUE, tolerance = 1E-3)

  #Hessian
  if (!requireNamespace("numDeriv", quietly = TRUE)){skip("Need numDeriv package")}
  numderiv <- numDeriv::hessian(ppill_r, u, beta0 = beta0, AL = AL, bL = bL)
  hess <- lltape$Hes(u, theta)
  dim(hess) <- c(p, p)
  expect_equal(hess, numderiv, ignore_attr = TRUE)

  #wrt u
  lltape_theta <- tape_swap(lltape)
  expect_equal(ppill_r(u, beta0, AL, bL), lltape_theta$eval(theta, u), ignore_attr = TRUE)

  #gradiant
  ppill_r_swap <- function(theta, u){
    pars <- ppi_parammats(theta)
    ppill_r(u, pars$beta, pars$AL, pars$bL)
  }
  numderiv <- numericDeriv(quote(ppill_r_swap(theta, u)), c("theta"))
  expect_equal(attr(numderiv,"gradient"), lltape_theta$Jac(theta, u),
               ignore_attr = TRUE, tolerance = 1E-7)

  #Hessian
  if (!requireNamespace("numDeriv", quietly = TRUE)){skip("Need numDeriv package")}
  numderiv <- numDeriv::hessian(ppill_r_swap, theta, u = u)
  hess <- lltape_theta$Hes(theta, u)
  dim(hess) <- c(length(theta), length(theta))
  expect_equal(hess, numderiv, ignore_attr = TRUE, tolerance = 1E-5)
})

# test Jacobian of ll function using numerical differentiation
test_that("ppi likelihood, Jacobian, Hessian for sphere matches numerical estimates", {
  psphere <- manifoldtransform("sim", "sqrt", "sph")
  lltape <- tapell("ppi", u, NA * theta, tran = psphere$tran, function(n){theta})

  # wrt u
  expect_equal(ppill_r_S(u, beta0, AL, bL), lltape$eval(u, theta), ignore_attr = TRUE)

  #gradiant
  numderiv <- numericDeriv(quote(ppill_r_S(u, beta0, AL, bL)), c("u"))
  expect_equal(attr(numderiv,"gradient"), lltape$Jac(u, theta), ignore_attr = TRUE, tolerance = 1E-3)

  #Hessian
  if (!requireNamespace("numDeriv", quietly = TRUE)){skip("Need numDeriv package")}
  numderiv <- numDeriv::hessian(ppill_r_S, u, beta0 = beta0, AL = AL, bL = bL)
  hess <- lltape$Hes(u, theta)
  dim(hess) <- c(p, p)
  expect_equal(hess, numderiv, ignore_attr = TRUE)

  #wrt theta
  lltape_theta <- tape_swap(lltape)
  expect_equal(ppill_r_S(u, beta0, AL, bL), lltape_theta$eval(theta, u), ignore_attr = TRUE)

  #gradiant
  ppill_r_S_swap <- function(theta, z){
    pars <- ppi_parammats(theta)
    ppill_r_S(z, pars$beta, pars$AL, pars$bL)
  }
  numderiv <- numericDeriv(quote(ppill_r_S_swap(theta, u)), c("theta"))
  expect_equal(attr(numderiv,"gradient"), lltape_theta$Jac(theta, u),
               ignore_attr = TRUE, tolerance = 1E-7)

  #Hessian
  if (!requireNamespace("numDeriv", quietly = TRUE)){skip("Need numDeriv package")}
  numderiv <- numDeriv::hessian(ppill_r_S_swap, theta, z = u)
  hess <- lltape_theta$Hes(theta, u)
  dim(hess) <- c(length(theta), length(theta))
  expect_equal(hess, numderiv, ignore_attr = TRUE, tolerance = 1E-5)
})


test_that("dirichlet ll evaluation and Jacobian matches expected", {
  beta = c(-0.3, -0.1, 3)
  set.seed(1234)
  u <- MCMCpack::rdirichlet(1, beta+1)

  dirichlet_r <- function(u, beta){sum(beta * log(u))}

  psimplex <- manifoldtransform("sim", "identity", "sim")
  lltape <- tapell("dirichlet", u, NA * beta, tran = psimplex$tran, function(n){beta})
  #forward0
  expect_equal(dirichlet_r(u, beta), lltape$eval(u, beta), ignore_attr = TRUE)

  #gradiant wrt u
  numderiv <- numericDeriv(quote(dirichlet_r(u, beta)), c("u"))
  expect_equal(attr(numderiv,"gradient"), lltape$Jac(u, beta), ignore_attr = TRUE, tolerance = 1E-5)

  #gradient wrt beta
  lltape_theta <- tape_swap(lltape)
  numderiv <- numericDeriv(quote(dirichlet_r(u, beta)), c("beta"))
  expect_equal(attr(numderiv,"gradient"), lltape_theta$Jac(beta, u), ignore_attr = TRUE, tolerance = 1E-5)
})
