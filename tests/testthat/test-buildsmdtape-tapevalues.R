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
  expect_equal(tapes1$uld_reembed$eval(ueval, thetaeval), tapes2$uld_reembed$eval(ueval, thetaeval))
  expect_equal(tapes1$uld_reembed$Jac(ueval, thetaeval),  tapes2$uld_reembed$Jac(ueval, thetaeval))
  expect_equal(tapes1$uld_reembed$Hes(ueval, thetaeval),  tapes2$uld_reembed$Hes(ueval, thetaeval))

  expect_equal(tapes1$smi$eval(thetaeval, ueval), tapes2$smi$eval(thetaeval, ueval))
  expect_equal(tapes1$smi$Jac(thetaeval, ueval),  tapes2$smi$Jac(thetaeval, ueval))
  expect_equal(tapes1$smi$Hes(thetaeval, ueval),  tapes2$smi$Hes(thetaeval, ueval))
  return(NULL)
  }

  #Sphere with minsq
  tapes1 <- tape_smi(manifold = "sph",
               uld = tape_uld_inbuilt("dirichlet", x = u1),
               transform = "sqrt",
               bdryw = tape_bdryw_inbuilt("minsq", u1, acut = acut)
               )
  tapes2 <- tape_smi(manifold = "sph",
                     uld = tape_uld_inbuilt("dirichlet", x = u2),
                     transform = "sqrt", 
                     bdryw = tape_bdryw_inbuilt("minsq", u2, acut = acut))

  compare(ueval, thetaeval, tapes1, tapes2)

  #Sphere, prodsq
  tapes1 <- tape_smi(manifold = "sph",
               uld = tape_uld_inbuilt("dirichlet", x = u1),
               transform = "sqrt",
               bdryw = tape_bdryw_inbuilt("prodsq", u1, acut = acut)
               )
  tapes2 <- tape_smi(manifold = "sph",
               uld = tape_uld_inbuilt("dirichlet", x = u2),
               transform = "sqrt",
               bdryw = tape_bdryw_inbuilt("prodsq", u2, acut = acut)
               )
  compare(ueval, thetaeval, tapes1, tapes2)

  #Simplex
  tapes1 <- tape_smi(manifold = "sim",
               uld = tape_uld_inbuilt("dirichlet", x = u1),
               bdryw = tape_bdryw_inbuilt("minsq", u1, acut = acut)
               )
  tapes2 <- tape_smi(manifold = "sim",
               uld = tape_uld_inbuilt("dirichlet", x = u2),
               bdryw = tape_bdryw_inbuilt("minsq", u2, acut = acut)
               )
  compare(ueval, thetaeval, tapes1, tapes2)

  #Simplex, prodsq
  tapes1 <- tape_smi(manifold = "sim",
               uld = tape_uld_inbuilt("dirichlet", x = u1),
               bdryw = tape_bdryw_inbuilt("prodsq", u1, acut = acut)
               )
  tapes2 <- tape_smi(manifold = "sim",
               uld = tape_uld_inbuilt("dirichlet", x = u2),
               bdryw = tape_bdryw_inbuilt("prodsq", u2, acut = acut)
               )
  compare(ueval, thetaeval, tapes1, tapes2)
})

