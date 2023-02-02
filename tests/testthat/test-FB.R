test_that("Fisher-Bingham likelihood runs and matches R code", {
  set.seed(321)
  p <- 3
  A <- matrix(NA, ncol = p, nrow = p)
  A[upper.tri(A)] <- runif(sum(upper.tri(A)))
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- c(runif(p-1), NA)
  diag(A)[2] <- 0
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
  m <- runif(p)
  m <- m / sqrt(sum(m^2))
  k <- 2
  theta <- FB_mats2theta(k, m, A)
  expect_equal(FB_theta2mats(theta),
    list(k = k,
         m = m,
         km = k*m,
         A = A))


  set.seed(123)
  sample <- matrix(runif(2 * p, -10, 10), nrow = 2)
  sample <- sample / sqrt(rowSums(sample^2))
  stopifnot(all(abs(sqrt(rowSums(sample^2)) - 1) < 1E-5))

  pman <- manifoldtransform("Snative")
  lltape <- ptapell(sample[1,], seq.int(1, length.out = length(theta)), llname = "FB", pman,
                    fixedtheta = rep(FALSE, length(theta)), verbose = FALSE)

  expect_equal(pForward0(lltape, sample[1, ], theta), log(qdFB(sample[1, ], k, m, A)),
               ignore_attr = TRUE)## very important to check a tape
  #deriv wrt u
  u <- sample[1, ]
  expect_equal(pJacobian(lltape, sample[1, ], theta),
               attr(numericDeriv(quote(log(qdFB(u, k, m, A))), "u"), "gradient"),
               ignore_attr = TRUE, tolerance = 1E-5)
  #deriv wrt theta
  lltape_theta <- swapDynamic(lltape, theta, sample[1, ])
  qdFB_theta <- function(u, theta){
    mats <- FB_theta2mats(theta)
    qdFB(u, mats$k, mats$m, mats$A)
  }
  expect_equal(pJacobian(lltape_theta, theta, u),
    attr(numericDeriv(quote(log(qdFB_theta(u, theta))), "theta"), "gradient"),
    ignore_attr = TRUE, tolerance = 1E-5)
})

test_that("rfb() simulation for diagonal matricies via Bingham() fitting", {
  skip_on_cran() #rfb() tested when fitting tested
  p <- 3
  A <- diag(c(30, 1, -31))
  sample <- Directional::rfb(1000, 1E-10, c(1, 0, 0), -A)
  est <- Bingham_full(sample)

  v010tonth <- rotationmatrix(c(1, rep(0, p-1)), c(0, 1, 0))
  rotatedA <- v010tonth %*% A %*% solve(v010tonth) #rotate A so that size of diagonal is increasing
  expect_lt_v(abs(est$A - rotatedA)[-(p*p)], 3 * est$A_SE[-(p*p)])  #the index removal of p*p removes that final element of the diagonal

  sample <- Directional::rfb(1000, 1E-10, c(1, 0, 0), -rotatedA)
  est <- Bingham_full(sample, control = list(tol = 1E-15))
  expect_lt_v(abs(est$A - A)[-(p*p)], 3 * est$A_SE[-(p*p)])  #the index removal of p*p removes that final element of the diagonal

  # try out FB estimation
  pman <- manifoldtransform("Snative")
  thetaFB <- FB_mats2theta(1E-10, c(1, 0, 0), A)
  tapes <- buildsmotape("Snative",
                        "FB",
                        utape = sample[1, ],
                        usertheta = rep(NA, length(thetaFB)))

  est <- cppad_closed(tapes$smotape, Y=sample)

  expect_absdiff_lte_v(est$est, thetaFB, 3 * est$SE)
})

test_that("rfb() simulation for general symmetric matrices via Bingham() fitting", {
  skip_on_cran() #rfb() tested when fitting tested
  p <- 3
  set.seed(345)
  theta <- runif(p-1 + (p - 1) * p/2)
  A <- Bingham_theta2Amat(theta)
  # diag(A) <- c(30, 1, NA) * diag(A)
  diag(A)[p] <- -sum(diag(A)[-p], na.rm = TRUE)
  theta <- Bingham_Amat2theta(A)
  A_es <- eigen(A)

  #simulate samples with a very small weighting on the concentration for von-mises
  set.seed(123)
  m <- c(1, 0, 0)
  B <- Directional::rotation(c(0, 1, 0), m)
  sample <- Directional::rfb(1000, 1E-8, m, -solve(B) %*% A %*% B)
  est <- Bingham_full(sample)
  estA_es <- eigen(est$A)

  #results have the same eigen values
  expect_equal(estA_es$values, A_es$values, tolerance = 0.5)
  #matrix values match
  expect_lt_v(abs(est$A - A)[-(p*p)], 3 * est$A_SE[-(p*p)])

  # note that parameter estimates have suprisingly large margin given that there is 1000 samples

  # try out FB estimation
  thetaFB <- FB_mats2theta(1E-8, m, A)
  tapes <- buildsmotape("Snative",
                        "FB",
                        utape = sample[1, ],
                        usertheta = rep(NA, length(thetaFB)))
                        
  est <- cppad_closed(tapes$smotape, Y=sample)

  expect_absdiff_lte_v(est$est, thetaFB, 3 * est$SE)
})


test_that("rfb() simulation for general symmetric matrices fitting", {
  p <- 3
  set.seed(111)
  theta <- runif(p-1 + (p - 1) * p/2 + p, -10, 10)
  thetamats <- FB_theta2mats(theta)
  # lapply(thetamats, round, 2) #nice view

  #simulate
  set.seed(12345)
  B <- Directional::rotation(c(0, 1, 0), thetamats$m)
  sample <- Directional::rfb(100000, thetamats$k, thetamats$m, -solve(B) %*% thetamats$A %*% B)

  #estimate
  tapes <- buildsmotape("Snative",
                        "FB",
                        utape = sample[1, ],
                        usertheta = rep(NA, length(theta)))
                        
  est <- cppad_closed(tapes$smotape, Y=sample)

  expect_absdiff_lte_v(est$est, theta, 3 * est$SE)
  # whooo it is working! But sample huge and SEs are still huge - it is almost like the model is misspecified
  # lapply(FB_theta2mats(est$par), round, 2)
  # lapply(thetamats, round, 2)
})


test_that("FB() fits for p = 3", {
  skip_on_cran() #test with various fixed elements is sufficient
  p <- 3
  set.seed(111)
  theta <- runif(p-1 + (p - 1) * p/2 + p, -10, 10)
  thetamats <- FB_theta2mats(theta)

  #simulate
  set.seed(12345)
  sample <- simdd::rFisherBingham(10000, mu = thetamats$k * thetamats$m, Aplus = thetamats$A)

  #Fit
  estobj <- FB(sample, control = list(tol = 1E-10))
  expect_lt_v(abs(estobj$est$paramvec - theta), 3 * estobj$SE$paramvec)
})

test_that("FB() fits with various fixed elements", {
  p <- 3
  set.seed(111)
  theta <- runif(p-1 + (p - 1) * p/2 + p, -10, 10)
  thetamats <- FB_theta2mats(theta)

  #simulate
  set.seed(12345)
  sample <- simdd::rFisherBingham(10000, mu = thetamats$k * thetamats$m, Aplus = thetamats$A)

  #a fixed A element
  inA <- matrix(NA, nrow = p, ncol = p)
  inA[p, 1] <- inA[1, p] <- thetamats$A[1, p]
  estobj <- FB(sample, A = inA, control = list(tol = 1E-12))
  expect_lte_v(abs(estobj$est$A - thetamats$A)[-(p*p)], 3 * estobj$SE$A[-(p*p)])
  expect_equal(estobj$est$A[!is.na(inA)], thetamats$A[!is.na(inA)])
  expect_equal(estobj$SE$A[!is.na(inA)], 0 * thetamats$A[!is.na(inA)])

  #fixed final diagonal element should error
  inA <- matrix(NA, nrow = p, ncol = p)
  inA[p, p] <- thetamats$A[p, p]
  expect_error(FB(sample, A = inA, control = list(tol = 1E-15)))

  # a fixed Fisher element
  inkm <- rep(NA, p)
  inkm[p] <- thetamats$m[p] * thetamats$k
  estobj <- FB(sample, km = inkm, control = list(tol = 1E-12))
  expect_lte_v(abs(estobj$est$km - thetamats$km), 3 * estobj$SE$km + 1E-10)
  expect_equal(estobj$est$km[!is.na(inkm)], thetamats$km[!is.na(inkm)])
  expect_equal(estobj$SE$km[!is.na(inkm)], 0 * thetamats$km[!is.na(inkm)])
})

test_that("FB() with many fixed elements leads to smaller smobjgrad", {
  skip("Fixing many of the elements doesn't improve smobjgrad")
  p <- 3
  set.seed(111)
  theta <- runif(p-1 + (p - 1) * p/2 + p, -10, 10)
  thetamats <- FB_theta2mats(theta)

  #simulate
  set.seed(12345)
  sample <- simdd::rFisherBingham(10000, mu = thetamats$k * thetamats$m, Aplus = thetamats$A)

  #many fixed elements
  intheta <- theta
  intheta[8] <- NA
  tapes <- buildsmotape("Snative", "FB",
                        sample[1, ], intheta)
  smvals <- tape_smvalues_wsum(tapes$smotape, sample, theta[is.na(intheta)])
  sum(smvals$grad^2)
})
