test_that("Dirichlet with smo values and derivatives independent of tape", {
  u1 = matrix(c(0.001, 0.011, 1 - 0.01 - 0.011), nrow = 1)
  u2 = matrix(c(1,1,1), nrow = 1)
  theta1 = rep(-0.5, 3)
  theta2 = rep(0, 3)
  fixedtheta = rep(FALSE, 3)
  acut = 0.1

  psphere <- manifoldtransform("sqrt", "sph")

  #Sphere
  pdir1 <- ptapell(u1, theta1, "dirichlet", psphere$tran, fixedtheta = fixedtheta, verbose = FALSE)
  pdir2 <- ptapell(u2, theta2, "dirichlet", psphere$tran, fixedtheta = fixedtheta, verbose = FALSE)
  pdir1smo <- ptapesmo(u1, theta1, pdir1, psphere$tran, psphere$man, "minsq", acut = acut, verbose = FALSE)
  pdir2smo <- ptapesmo(u2, theta2, pdir2, psphere$tran, psphere$man, "minsq", acut = acut, verbose = FALSE)

  ueval <- matrix(c(0.4, 0.011, 1 - 0.4 - 0.011), nrow = 1)
  thetaeval <- c(-0.1, -0.5, 2)
  expect_equal(pForward0(pdir1, ueval, thetaeval), pForward0(pdir2, ueval, thetaeval))
  expect_equal(pJacobian(pdir1, ueval, thetaeval), pJacobian(pdir2, ueval, thetaeval))
  expect_equal(pHessian(pdir1, ueval, thetaeval), pHessian(pdir2, ueval, thetaeval))

  expect_equal(pForward0(pdir1smo, thetaeval, ueval), pForward0(pdir2smo, thetaeval, ueval))
  expect_equal(pJacobian(pdir1smo, thetaeval, ueval), pJacobian(pdir2smo, thetaeval, ueval))
  expect_equal(pHessian(pdir1smo, thetaeval, ueval), pHessian(pdir2smo, thetaeval, ueval))

  #Sphere, prodsq
  pdir1 <- ptapell(u1, theta1, "dirichlet", psphere, fixedtheta = fixedtheta, verbose = FALSE)
  pdir2 <- ptapell(u2, theta2, "dirichlet", psphere, fixedtheta = fixedtheta, verbose = FALSE)
  pdir1smo <- ptapesmo(u1, theta1, pdir1, psphere, "prodsq", acut = acut, verbose = FALSE)
  pdir2smo <- ptapesmo(u2, theta2, pdir2, psphere, "prodsq", acut = acut, verbose = FALSE)

  ueval <- matrix(c(0.4, 0.011, 1 - 0.4 - 0.011), nrow = 1)
  thetaeval <- c(-0.1, -0.5, 2)
  expect_equal(pForward0(pdir1, ueval, thetaeval), pForward0(pdir2, ueval, thetaeval))
  expect_equal(pJacobian(pdir1, ueval, thetaeval), pJacobian(pdir2, ueval, thetaeval))
  expect_equal(pHessian(pdir1, ueval, thetaeval), pHessian(pdir2, ueval, thetaeval))

  expect_equal(pForward0(pdir1smo, thetaeval, ueval), pForward0(pdir2smo, thetaeval, ueval))
  expect_equal(pJacobian(pdir1smo, thetaeval, ueval), pJacobian(pdir2smo, thetaeval, ueval))
  expect_equal(pHessian(pdir1smo, thetaeval, ueval), pHessian(pdir2smo, thetaeval, ueval))

  #Simplex
  psimplex <- manifoldtransform("identity", "sim")
  pdir1 <- ptapell(u1, theta1, "dirichlet", psimplex$tran, fixedtheta = fixedtheta, verbose = FALSE)
  pdir2 <- ptapell(u2, theta2, "dirichlet", psimplex$tran, fixedtheta = fixedtheta, verbose = FALSE)
  pdir1smo <- ptapesmo(u1, theta1, pdir1, psimplex$tran, psimplex$man, "minsq", acut = acut, verbose = FALSE)
  pdir2smo <- ptapesmo(u2, theta2, pdir2, psimplex$tran, psimplex$man, "minsq", acut = acut, verbose = FALSE)

  expect_equal(pForward0(pdir1, ueval, thetaeval), pForward0(pdir2, ueval, thetaeval))
  expect_equal(pJacobian(pdir1, ueval, thetaeval), pJacobian(pdir2, ueval, thetaeval))
  expect_equal(pHessian(pdir1, ueval, thetaeval), pHessian(pdir2, ueval, thetaeval))

  expect_equal(pForward0(pdir1smo, thetaeval, ueval), pForward0(pdir2smo, thetaeval, ueval))
  expect_equal(pJacobian(pdir1smo, thetaeval, ueval), pJacobian(pdir2smo, thetaeval, ueval))
  expect_equal(pHessian(pdir1smo, thetaeval, ueval), pHessian(pdir2smo, thetaeval, ueval))

  #Simplex, prodsq
  pdir1 <- ptapell(u1, theta1, "dirichlet", psimplex, fixedtheta = fixedtheta, verbose = FALSE)
  pdir2 <- ptapell(u2, theta2, "dirichlet", psimplex, fixedtheta = fixedtheta, verbose = FALSE)
  pdir1smo <- ptapesmo(u1, theta1, pdir1, psimplex, "prodsq", acut = acut, verbose = FALSE)
  pdir2smo <- ptapesmo(u2, theta2, pdir2, psimplex, "prodsq", acut = acut, verbose = FALSE)

  expect_equal(pForward0(pdir1, ueval, thetaeval), pForward0(pdir2, ueval, thetaeval))
  expect_equal(pJacobian(pdir1, ueval, thetaeval), pJacobian(pdir2, ueval, thetaeval))
  expect_equal(pHessian(pdir1, ueval, thetaeval), pHessian(pdir2, ueval, thetaeval))

  expect_equal(pForward0(pdir1smo, thetaeval, ueval), pForward0(pdir2smo, thetaeval, ueval))
  expect_equal(pJacobian(pdir1smo, thetaeval, ueval), pJacobian(pdir2smo, thetaeval, ueval))
  expect_equal(pHessian(pdir1smo, thetaeval, ueval), pHessian(pdir2smo, thetaeval, ueval))
})

