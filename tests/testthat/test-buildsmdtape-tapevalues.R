test_that("Dirichlet with smd values and derivatives independent of tape", {
  u1 = matrix(c(0.001, 0.011, 1 - 0.01 - 0.011), nrow = 1)
  u2 = matrix(c(1,1,1), nrow = 1)
  theta1 = rep(-0.5, 3)
  theta2 = rep(0, 3)
  fixedtheta = rep(FALSE, 3)
  acut = 0.1
  ueval <- matrix(c(0.4, 0.011, 1 - 0.4 - 0.011), nrow = 1)
  thetaeval <- c(-0.1, -0.5, 2)

  compare <- function(uval, thetaeval, tapes1, tapes2){
  expect_equal(tapes1$lltape$eval(ueval, thetaeval), tapes2$lltape$eval(ueval, thetaeval))
  expect_equal(tapes1$lltape$Jac(ueval, thetaeval), tapes2$lltape$Jac(ueval, thetaeval))
  expect_equal(tapes1$lltape$Hes(ueval, thetaeval), tapes2$lltape$Hes(ueval, thetaeval))

  expect_equal(tapes1$smdtape$eval(thetaeval, ueval), tapes2$smdtape$eval(thetaeval, ueval))
  expect_equal(tapes1$smdtape$Jac(thetaeval, ueval), tapes2$smdtape$Jac(thetaeval, ueval))
  expect_equal(tapes1$smdtape$Hes(thetaeval, ueval), tapes2$smdtape$Hes(thetaeval, ueval))
  return(NULL)
  }

  #Sphere with minsq
  tapes1 <- tape_smd("sim", "sqrt", "sph",
               "dirichlet", u1, rep(NA, 3), thetatape_creator = function(n){theta1},
               bdryw = "minsq", acut = acut)
  tapes2 <- tape_smd("sim", "sqrt", "sph",
               "dirichlet", u2, rep(NA, 3), thetatape_creator = function(n){theta2},
               bdryw = "minsq", acut = acut)

  compare(ueval, thetaeval, tapes1, tapes2)

  #Sphere, prodsq
  tapes1 <- tape_smd("sim", "sqrt", "sph",
               "dirichlet", u1, rep(NA, 3), thetatape_creator = function(n){theta1},
               bdryw = "prodsq", acut = acut)
  tapes2 <- tape_smd("sim", "sqrt", "sph",
               "dirichlet", u2, rep(NA, 3), thetatape_creator = function(n){theta2},
               bdryw = "prodsq", acut = acut)
  compare(ueval, thetaeval, tapes1, tapes2)

  #Simplex
  tapes1 <- tape_smd("sim", "identity", "sim",
               "dirichlet", u1, rep(NA, 3), thetatape_creator = function(n){theta1},
               bdryw = "minsq", acut = acut)
  tapes2 <- tape_smd("sim", "identity", "sim",
               "dirichlet", u2, rep(NA, 3), thetatape_creator = function(n){theta2},
               bdryw = "minsq", acut = acut)
  compare(ueval, thetaeval, tapes1, tapes2)

  #Simplex, prodsq
  tapes1 <- tape_smd("sim", "identity", "sim",
               "dirichlet", u1, rep(NA, 3), thetatape_creator = function(n){theta1},
               bdryw = "prodsq", acut = acut)
  tapes2 <- tape_smd("sim", "identity", "sim",
               "dirichlet", u2, rep(NA, 3), thetatape_creator = function(n){theta2},
               bdryw = "prodsq", acut = acut)
  compare(ueval, thetaeval, tapes1, tapes2)
})

