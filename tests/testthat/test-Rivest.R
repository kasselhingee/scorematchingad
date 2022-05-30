test_that("Rivest likelihood runs and matches R code", {
  set.seed(321)
  p <- 3
  A <- rsymmetricmatrix(p, -10, 10)
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
  A_es <- eigen(A)
  idx <- 1 #smallest in abs size
  k <- 2

  theta <- c(cdabyppi:::Bingham_Amat2theta(A), k, idx)

  set.seed(123)
  sample <- matrix(runif(2 * p, -10, 10), nrow = 2)
  sample <- sample / sqrt(rowSums(sample^2))
  stopifnot(all(abs(sqrt(rowSums(sample^2)) - 1) < 1E-5))

  pman <- pmanifold("Snative")
  lltape <- ptapell(sample[1,], theta, llname = "Rivest", pman,
                    fixedtheta = rep(FALSE, length(theta)), verbose = FALSE)

  u <- sample[2, ]
  expect_equal(pForward0(lltape, u, theta), log(cdabyppi:::qdRivest(u, k, A, idx)),
               ignore_attr = TRUE)
  # different eigenvectors
  expect_equal(pForward0(lltape, u, c(theta[-length(theta)], 2)), log(cdabyppi:::qdRivest(u, k, A, 2)),
               ignore_attr = TRUE)
  expect_equal(pForward0(lltape, u, c(theta[-length(theta)], 3)), log(cdabyppi:::qdRivest(u, k, A, 3)),
               ignore_attr = TRUE)

  #derive wrt u
  expect_equal(pJacobian(lltape, u, theta), cdabyppi:::lldRivest_du(u, k, A, idx),
               ignore_attr = TRUE)

  #deriv wrt theta
  llRivest <- function(theta){
    A <- cdabyppi:::Bingham_theta2Amat(theta[seq.int(1, length.out = p - 1 + (p-1)*p/2)])
    k <- theta[1 + p - 1 + (p-1)*p/2]
    idx <- theta[2 + p - 1 + (p-1)*p/2]
    out <- log(cdabyppi:::qdRivest(u, k, A, idx))
    return(out)
  }

  # test changing A gives different eigenvector - works as of 26 May, 2022, purely using PrintFor statement so can't test on it
  thetatest <- theta
  thetatest[c(1,3,4)] <- thetatest[c(1,3,4)] / 100
  thetatest[length(thetatest)] <- 3 #choose the largest rather than the smallest
  expect_equal(pForward0(lltape, u, thetatest), llRivest(thetatest),
               ignore_attr = TRUE)

  # test deriv wrt theta
  Rgradt <- numericDeriv(quote(llRivest(thetatest)), c("thetatest"))
  lltape_t <- swapDynamic(lltape, theta, sample[1, ])
  expect_equal(pJacobian(lltape_t, thetatest, u), attr(Rgradt, "gradient"),
               tolerance = 1E-5, ignore_attr = TRUE)
})

test_that("ll_Rivest modifies m to have the correct sign", {
  # the eigenvector with the second highest eigenvalue comes out of eigendecomposition with positive first element
  # the first element has been negative in the past
  p <- 3
  k <- -3.2
  evidx <- 2
  theta <- c(9.136667, -0.9333169, 7.8483809, 1.0287003, 0.8677053, k, evidx)
  mats <- Rivest_theta2mats(theta)
  # [,1]       [,2]       [,3]
  # [1,] 9.136667  7.8483809  1.0287003
  # [2,] 7.848381 -0.9333169 -0.8677053
  # [3,] 1.028700 -0.8677053 -8.2033500

  set.seed(123)
  sample <- matrix(runif(2 * p, -10, 10), nrow = 2)
  sample <- sample / sqrt(rowSums(sample^2))
  stopifnot(all(abs(sqrt(rowSums(sample^2)) - 1) < 1E-5))

  ltheta <- length(theta)
  thetafortape <- c(seq.int(1, length.out = ltheta-1), evidx)
  pman <- pmanifold("Snative")
  lltape <- ptapell(sample[1,], thetafortape, llname = "Rivest", pman,
                    fixedtheta = rep(FALSE, ltheta), verbose = FALSE)
  # for this set of theta, the eigen value has a positive first element!
  expect_equal( pForward0(lltape, sample[2, ], theta), log(qdRivest(sample[2, ], mats$k, mats$A, mats$evidx)),
               ignore_attr = TRUE)

  #derive wrt u
  expect_equal(pJacobian(lltape, sample[2, ], theta), cdabyppi:::lldRivest_du(sample[2, ], mats$k, mats$A, mats$evidx),
               ignore_attr = TRUE)

  #deriv wrt theta
  llRivest <- function(theta){
    A <- cdabyppi:::Bingham_theta2Amat(theta[seq.int(1, length.out = p - 1 + (p-1)*p/2)])
    k <- theta[1 + p - 1 + (p-1)*p/2]
    idx <- theta[2 + p - 1 + (p-1)*p/2]
    out <- log(cdabyppi:::qdRivest(sample[2, ], k, A, idx))
    return(out)
  }

  Rgradt <- numericDeriv(quote(llRivest(theta)), c("theta"))
  lltape_t <- swapDynamic(lltape, theta, sample[1, ])
  expect_equal(pJacobian(lltape_t, theta, sample[2, ]), attr(Rgradt, "gradient"),
               tolerance = 1E-5, ignore_attr = TRUE)
})

test_that("taped ll_Rivest is the same for differrent tape values", {
  set.seed(321)
  p <- 3
  A <- rsymmetricmatrix(p, -10, 10)
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
  A_es <- eigen(A)
  evidx <- 1 #smallest in abs size
  k <- 2

  theta <- c(cdabyppi:::Bingham_Amat2theta(A), k, evidx)

  set.seed(123)
  sample <- matrix(runif(2 * p, -10, 10), nrow = 2)
  sample <- sample / sqrt(rowSums(sample^2))
  stopifnot(all(abs(sqrt(rowSums(sample^2)) - 1) < 1E-5))

  ltheta <- length(theta)
  # fix both k and evidx correctly
  thetafortape <- c(seq.int(1, length.out = ltheta-1), theta[ltheta])
  thetafortape[ltheta - 1] <- theta[ltheta - 1]
  fixedtheta <- rep(FALSE, ltheta)

  pman <- pmanifold("Snative")
  lltape_A <- ptapell(sample[1,], thetafortape, llname = "Rivest", pman,
                    fixedtheta = fixedtheta, verbose = FALSE)
  expect_equal(pForward0(lltape_A, sample[2, ], theta), log(qdRivest(sample[2, ], k, A, evidx)),
               ignore_attr = TRUE) #FAILS

  lltape <- ptapell(sample[1,], theta, llname = "Rivest", pman,
                    fixedtheta = fixedtheta, verbose = FALSE)
  expect_equal(pForward0(lltape, sample[2, ], theta), log(qdRivest(sample[2, ], k, A, evidx)),
               ignore_attr = TRUE)
  expect_equal(pForward0(lltape, sample[1, ], theta), log(qdRivest(sample[1, ], k, A, evidx)),
               ignore_attr = TRUE)

  expect_equal(pForward0(lltape, sample[2, ], theta), pForward0(lltape_A, sample[2, ], theta))
  expect_equal(pForward0(lltape, sample[2, ], thetafortape), pForward0(lltape_A, sample[2, ], thetafortape))

}) #failing

test_that("Taped Rivest smo has derivatives consistent with numerical differentation", {
  set.seed(321)
  p <- 3
  A <- rsymmetricmatrix(p, -10, 10)
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
  A_es <- eigen(A)
  evidx <- 1 #smallest in abs size
  k <- 2

  theta <- c(cdabyppi:::Bingham_Amat2theta(A), k, evidx)

  set.seed(123)
  sample <- matrix(runif(2 * p, -10, 10), nrow = 2)
  sample <- sample / sqrt(rowSums(sample^2))
  stopifnot(all(abs(sqrt(rowSums(sample^2)) - 1) < 1E-5))

  ltheta <- length(theta)
  # fix both k and evidx correctly
  thetafortape <- c(seq.int(1, length.out = ltheta-1), theta[ltheta])
  thetafortape[ltheta - 1] <- theta[ltheta - 1]
  fixedtheta <- rep(FALSE, ltheta)

  pman <- pmanifold("Snative")
  lltape <- ptapell(sample[1,], thetafortape, llname = "Rivest", pman,
                    fixedtheta = fixedtheta, verbose = FALSE)
  smotape <- ptapesmo(sample[1,], thetafortape[!fixedtheta],
                      lltape, pman, "ones", 1, verbose = FALSE)

  val <- pForward0(smotape, theta, sample[2, ])
  rgrad <- numericDeriv(quote(pForward0(smotape, theta, sample[2, ])), "theta")
  stopifnot(abs(rgrad - val) < 1E-5)

  expect_equal(pJacobian(smotape, theta, sample[2, ]), attr(rgrad, "gradient"),
               ignore_attr = TRUE, tolerance = 1E-5)
})

test_that("Rivest() fits correctly with evidx given correctly", {
  p <- 3
  set.seed(1245)
  A <- rsymmetricmatrix(p, -10, 10)
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint for Bingham
  k <- -3.2
  evidx <- 2
  sample <- rRivest(1000, k, A, evidx)
  theta <- Rivest_mats2theta(k, A, evidx)
  control <- list(tol = 1E-10)

  ltheta <- p-1 + (p - 1) * p/2 + 2
  expect_equal(ltheta, length(theta))
  # fix both k and evidx correctly
  thetafortape <- c(seq.int(1, length.out = ltheta-1), theta[ltheta])
  thetafortape[ltheta - 1] <- theta[ltheta - 1]
  fixedtheta <- c(rep(FALSE, ltheta-2), FALSE, TRUE)

  pman <- pmanifold("Snative")
  lltape <- ptapell(sample[1,], thetafortape, llname = "Rivest", pman,
                    fixedtheta = fixedtheta, verbose = FALSE)
  # for this set of theta, the eigen value has a positive first element!
  expect_equal(pForward0(lltape, sample[2, ], theta[!fixedtheta]), log(qdRivest(sample[2, ], k, A, evidx)),
               ignore_attr = TRUE)


  smotape <- ptapesmo(sample[1,], thetafortape[!fixedtheta],
                        lltape, pman, "ones", 1, verbose = FALSE)

  # first there is something strange that the gradient of smotape is very high at the correct value of the theta
  smobj(smotape, theta[!fixedtheta], sample)
  smobjgrad(smotape, theta[!fixedtheta], sample)
  pJacobian(smotape, theta[!fixedtheta], sample[1, ])

  smobjhess(smotape, theta[!fixedtheta], sample)
  smestSE(smotape, theta[!fixedtheta], sample)

  # minimise without using smobjgrad
  out <- Rcgmin::Rcgmin(par = thetafortape[!fixedtheta],
                        fn = function(theta){smobj(smotape, theta, sample)},
                        control = control)
  outSE <- smestSE(smotape, out$par, sample)
  expect_lt_v(abs(out$par - theta[!fixedtheta]), 3 * outSE)
})

test_that("Rivest() fits correctly when starting values", {
  p <- 3
  set.seed(1245)
  A <- rsymmetricmatrix(p, -10, 10)
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint for Bingham
  k <- 1E-5
  evidx <- 2
  sample <- rRivest(1000, k, A, evidx)
  theta <- Rivest_mats2theta(k, A, evidx)
  control <- list(tol = 1E-10)

  ltheta <- p-1 + (p - 1) * p/2 + 2
  expect_equal(ltheta, length(theta))
  # fix both k and evidx correctly
  thetafortape <- c(seq.int(1, length.out = ltheta-1), theta[ltheta])
  thetafortape[ltheta - 1] <- theta[ltheta - 1]
  fixedtheta <- c(rep(FALSE, ltheta-2), FALSE, TRUE)

  pman <- pmanifold("Snative")
  lltape <- ptapell(sample[1,], thetafortape, llname = "Rivest", pman,
                    fixedtheta = fixedtheta, verbose = FALSE)
  # for this set of theta, the eigen value has a positive first element!
  expect_equal(pForward0(lltape, sample[2, ], theta[!fixedtheta]), log(qdRivest(sample[2, ], k, A, evidx)),
               ignore_attr = TRUE)


  smotape <- ptapesmo(sample[1,], thetafortape[!fixedtheta],
                      lltape, pman, "ones", 1, verbose = FALSE)

  # first there is something strange that the gradient of smotape is very high at the correct value of the theta
  smobj(smotape, theta[!fixedtheta], sample)
  smobjgrad(smotape, theta[!fixedtheta], sample)
  pJacobian(smotape, theta[!fixedtheta], sample[1, ])

  smobjhess(smotape, theta[!fixedtheta], sample)
  smestSE(smotape, theta[!fixedtheta], sample)

  # minimise without using smobjgrad
  out <- Rcgmin::Rcgmin(par = thetafortape[!fixedtheta],
                        fn = function(theta){smobj(smotape, theta, sample)},
                        control = control)
  outSE <- smestSE(smotape, out$par, sample)
  expect_lt_v(abs(out$par - theta[!fixedtheta]), 3 * outSE)
})



  sminfo <- smest(smotape, thetafortape[-ltheta], sample,
                  control = control)
  thetamat <- FB_theta2mats(sminfo$par)
  smestSE(smotape, sminfo$par, sample)

})
