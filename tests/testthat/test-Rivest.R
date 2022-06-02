test_that("Rivest likelihood and derives match R code when taped using same theta", {
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

  #double derive wrt u
  expect_equal(pHessian(lltape, u, theta), as.vector(cdabyppi:::lldRivest_dudu(u, k, A, idx)),
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

  # test deriv wrt theta at theta
  Rgradt <- numericDeriv(quote(llRivest(theta)), c("theta"))
  lltape_t <- swapDynamic(lltape, theta, sample[1, ])
  expect_equal(pJacobian(lltape_t, theta, u), attr(Rgradt, "gradient"),
               tolerance = 1E-5, ignore_attr = TRUE)

  # test deriv wrt theta at thetatest
  Rgradt <- numericDeriv(quote(llRivest(thetatest)), c("thetatest"))
  lltape_t <- swapDynamic(lltape, theta, sample[1, ])
  expect_equal(pJacobian(lltape_t, thetatest, u), attr(Rgradt, "gradient"),
               tolerance = 1E-5, ignore_attr = TRUE)
})

test_that("Rivest_theta2FBtheta() leads to correct likelihood values", {
  set.seed(321)
  p <- 3
  A <- rsymmetricmatrix(p, -10, 10)
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
  A_es <- eigen(A)
  evidx <- 1 #smallest in abs size
  k <- -2

  theta <- c(cdabyppi:::Bingham_Amat2theta(A), k, evidx)

  set.seed(123)
  sample <- matrix(runif(2 * p, -10, 10), nrow = 2)
  sample <- sample / sqrt(rowSums(sample^2))

  FBtheta <- Rivest_theta2FBtheta(theta)
  FBmats <- FB_theta2mats(FBtheta)
  expect_equal(FBmats$A, A)
  expect_equal(FBmats$k, abs(k))
  expect_equal(qdFB(sample[1, ], FBmats$k, FBmats$m, FBmats$A), qdRivest(sample[1, ], k, A, evidx))
  expect_equal(qdFB(sample[2, ], FBmats$k, FBmats$m, FBmats$A), qdRivest(sample[2, ], k, A, evidx))
})

test_that("ll_Rivest taped on correct theta gives correct results for many measurments", {
  p <- 3
  set.seed(124567)
  A <- rsymmetricmatrix(p, -10, 10)
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint for Bingham
  k <- -3.2
  evidx <- 2
  sample <- rRivest(10000, k, A, evidx)
  theta <- Rivest_mats2theta(k, A, evidx)

  #check that tape at true theta gives good results
  tapes <- buildsmotape("Snative", "Rivest",
                        sample[1, ], theta * NA,
                        thetatape_creator = function(n){theta[1:n]})
  stopifnot(all(tapes$info$tapedtheta == theta))
  llcppad <- unlist(lapply(1:nrow(sample), function(i){pForward0(tapes$lltape, sample[i, ], theta)}))
  llR <- unlist(lapply(1:nrow(sample), function(i){log(qdRivest(sample[i, ], k, A, evidx))}))
  expect_true(all(abs(llcppad - llR) < 1E-10))

  # and gradient wrt to u!
  dllcppad <- do.call(rbind, lapply(1:nrow(sample), function(i){pJacobian(tapes$lltape, sample[i, ], theta)}))
  dllR <- do.call(rbind, lapply(1:nrow(sample), function(i){t(lldRivest_du(sample[i, ], k, A, evidx))}))
  expect_true(all(abs(dllcppad - dllR) < 1E-10))

  # and Hessian wrt to u
  ddllcppad <- do.call(rbind, lapply(1:nrow(sample), function(i){pHessian(tapes$lltape, sample[i, ], theta)}))
  ddllR <- do.call(rbind, lapply(1:nrow(sample), function(i){as.vector(lldRivest_dudu(sample[i, ], k, A, evidx))}))
  expect_true(all(abs(ddllcppad - ddllR) < 1E-10))
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

test_that("taped ll_Rivest results are sensitive to the tape values", {
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
  expect_error(expect_equal(pForward0(lltape_A, sample[2, ], theta), log(qdRivest(sample[2, ], k, A, evidx)),
               ignore_attr = TRUE))

  lltape <- ptapell(sample[1,], theta, llname = "Rivest", pman,
                    fixedtheta = fixedtheta, verbose = FALSE)
  expect_equal(pForward0(lltape, sample[2, ], theta), log(qdRivest(sample[2, ], k, A, evidx)),
               ignore_attr = TRUE)
  expect_equal(pForward0(lltape, sample[1, ], theta), log(qdRivest(sample[1, ], k, A, evidx)),
               ignore_attr = TRUE)

  expect_error(expect_equal(pForward0(lltape, sample[2, ], theta), pForward0(lltape_A, sample[2, ], theta)))
  expect_equal(pForward0(lltape, sample[2, ], thetafortape), pForward0(lltape_A, sample[2, ], thetafortape))

})

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

test_that("Rivest matches Bingham when Fisher concentration is zero", {
  p <- 3
  set.seed(1245)
  A <- rsymmetricmatrix(p, -10, 10)
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint for Bingham
  k <- 0
  evidx <- 2
  sample <- rRivest(10000, k, A, evidx)
  # sample <- rBingham(10000, A)

  # Rivest ll funs etc
  Rtheta <- Rivest_mats2theta(k, A, evidx)
  #check that tape at true theta gives good results
  Rtapes <- buildsmotape("Snative", "Rivest",
                        sample[1, ], Rtheta * NA,
                        thetatape_creator = function(n){Rtheta[1:n]})

  Btheta <- Bingham_Amat2theta(A)
  Btapes <- buildsmotape("Snative", "Rivest",
                                sample[1, ], Btheta * NA,
                                thetatape_creator = function(n){Btheta[1:n]})
  #preliminary is ll evaulation - should be correct as this is tested against R code elsewhere in this file
  expect_equal(pForward0(Btapes$lltape, sample[2, ], Btheta), pForward0(Rtapes$lltape, sample[2, ], Rtheta))

  expect_equal(pForward0(Btapes$smotape, Btheta, sample[1, ]), pForward0(Rtapes$smotape, Rtheta, sample[1, ]))
  Rsmo <- unlist(lapply(1:nrow(sample), function(i){pForward0(Rtapes$smotape, Rtheta, sample[i, ])}))
  Bsmo <- unlist(lapply(1:nrow(sample), function(i){pForward0(Btapes$smotape, Btheta, sample[i, ])}))
  expect_true(all(abs(Rsmo - Bsmo) < 1E-10))

  #expected the the smo grad wrt theta is different for some measurements, but it isn't!!!
  Rsmograd <- do.call(rbind, lapply(1:nrow(sample), function(i){pJacobian(Rtapes$smotape, Rtheta, sample[i, ])}))
  Bsmograd <- do.call(rbind, lapply(1:nrow(sample), function(i){pJacobian(Btapes$smotape, Btheta, sample[i, ])}))
  expect_true(all(abs(Rsmograd[, 1:5] - Bsmograd) < 1E-10))

  expect_equal(smobjgrad(Rtapes$smotape, Rtheta, sample)[1:5], smobjgrad(Btapes$smotape, Btheta, sample))
})

test_that("Rivest matches Fisher-Bingham on smoval, but not smograd", {
  p <- 3
  set.seed(1245)
  A <- rsymmetricmatrix(p, -10, 10)
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint for Bingham
  k <- 3.2
  evidx <- 2
  A_es <- eigen(A)
  evalorder <- order(abs(A_es$values), decreasing = FALSE)
  m <- A_es$vectors[, evalorder[evidx]]
  if (m[1] < 0){m <- -m}
  sample <- rRivest(10000, k, A, evidx)

  # Rivest ll funs etc
  Rtheta <- Rivest_mats2theta(k, A, evidx)
  #check that tape at true theta gives good results
  Rtapes <- buildsmotape("Snative", "Rivest",
                         sample[1, ], Rtheta * NA,
                         thetatape_creator = function(n){Rtheta[1:n]})

  FBtheta <- FB_mats2theta(k, sign(k) * m, A)
  FBtapes <- buildsmotape("Snative", "FB",
                         sample[1, ], FBtheta * NA,
                         thetatape_creator = function(n){FBtheta[1:n]})

  expect_equal(pForward0(FBtapes$smotape, FBtheta, sample[1, ]), pForward0(Rtapes$smotape, Rtheta, sample[1, ]))
  Rsmo <- unlist(lapply(1:nrow(sample), function(i){pForward0(Rtapes$smotape, Rtheta, sample[i, ])}))
  FBsmo <- unlist(lapply(1:nrow(sample), function(i){pForward0(FBtapes$smotape, FBtheta, sample[i, ])}))
  expect_true(all(abs(Rsmo - FBsmo) < 1E-10))

  #expected the the smo grad wrt theta is different because the m vector depends on other theta in the Rivest distribution
  Rsmograd <- do.call(rbind, lapply(1:nrow(sample), function(i){pJacobian(Rtapes$smotape, Rtheta, sample[i, ])}))
  FBsmograd <- do.call(rbind, lapply(1:nrow(sample), function(i){pJacobian(FBtapes$smotape, FBtheta, sample[i, ])}))
  expect_false(all(abs(Rsmograd[, 1:5] - FBsmograd[, 1:5]) < 1E-10))
})

test_that("Fitting via smo val and FB evaluation", {
  p <- 3
  set.seed(1245)
  A <- rsymmetricmatrix(p, -10, 10)
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint for Bingham
  k <- -3.2
  evidx <- 2
  sample <- rRivest(10000, k, A, evidx)
  Rtheta <- Rivest_mats2theta(k, A, evidx)

  #first check that FB gradient and value is low
  FBtheta <- Rivest_theta2FBtheta(Rtheta)
  FBtapes <- buildsmotape("Snative", "FB",
                        sample[1, ], FBtheta,
                        thetatape_creator = function(n){FBtheta[1:n]})

})

test_that("Fitting when each unique theta is newly taped works for fixed evidx", {
  p <- 3
  set.seed(1245)
  A <- rsymmetricmatrix(p, -10, 10)
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint for Bingham
  k <- 0
  evidx <- 2
  # sample <- rRivest(10000, k, A, evidx)
  sample <- rBingham(10000, A)
  theta <- Rivest_mats2theta(k, A, evidx)
  control <- list(tol = 1E-10)

  #check that tape at true theta gives good results
  tapes <- buildsmotape("Snative", "Rivest",
               sample[1, ], theta * NA,
               thetatape_creator = function(n){theta[1:n]})
  stopifnot(abs(pForward0(tapes$lltape, sample[1, ], theta) - log(qdRivest(sample[1, ], k, A, evidx))) < 1E-10)
  smobj(tapes$smotape, theta, sample)
  smobjgrad(tapes$smotape, theta, sample)

  smoRivest_pertheta_peru <- function(theta, ufortaping = sample[1, ]){
    pman <- pmanifold("Snative")
    lltape <- ptapell(ufortaping, theta, llname = "Rivest", pman,
                      fixedtheta = c(rep(TRUE, length(theta) -1), FALSE), verbose = FALSE)
    smotape <- ptapesmo(ufortaping, theta[length(theta)],
                        lltape, pman, "ones", 1, verbose = FALSE)
    function(u){pForward0(smotape, theta[length(theta)], u)}
  }

  smoRivest_pertheta <- function(theta, sample, ufortaping = sample[1, ]){
    ufun <- smoRivest_pertheta_peru(theta, ufortaping = ufortaping)
    sc_perpt <- lapply(1:nrow(sample), function(i){
      scobj <- ufun(sample[i,])
      return(scobj)
    })
    scmo <- mean(unlist(sc_perpt))
    return(scmo)
  }

  # minimise without using smobjgrad
  ltheta <- p-1 + (p - 1) * p/2 + 2
  starttheta <- c(seq.int(1, length.out = ltheta-1), theta[ltheta])
  out <- Rcgmin::Rcgmin(par = starttheta,
                        fn = function(theta){smoRivest_pertheta(theta, sample)},
                        control = control)
  expect_equal(out$par, theta, tolerance = 1E-1, ignore_attr = TRUE)
})
