test_that("von-Mises Fisher likelihood runs and fits", {
  set.seed(123)
  theta <-  3 * c(1, -1) / sqrt(2)
  sample <- movMF::rmovMF(1000, theta)
  thetatape <- c(1, -1)

  p <- length(theta)
  intheta <- c(NA, rep(0, p - 1))
  tapes <- buildsmotape("Snative", "vMF",
                        rep(1, p)/sqrt(p), rep(NA, p),
                        weightname = "ones",
                        verbose = FALSE)
  expect_equal(pForward0(tapes$lltape, sample[1, ], theta), sum(sample[1, ]  * theta)) ## very important to check a tape
  out <- smest(tapes$smotape, thetatape, sample, control = list(tol = 1E-10))
  cdabyppi:::expect_lt_v(abs(out$par - theta), 3 * out$SE)
})

test_that("vMF_Mardia() function works on data around north pole", {
  set.seed(123)
  k <- 3
  m <- c(1, 0)
  km <-  k * m
  warning("switch to simulating vMF using Directional")
  sample <- movMF::rmovMF(1000, km)
  out <- vMF(sample, method = "Mardia")
  expect_equal(out$m , m, tolerance = 1E-1) #moment estimate part
  cdabyppi:::expect_lt_v(abs(out$k - k), 3 * out$SE$k)
})

test_that("vMF_Mardia() function works for data centred elsewhere", {
  set.seed(123)
  k <- 3
  m <- c(1, -1) / sqrt(2)
  km <-  k * m
  sample <- movMF::rmovMF(1000, km)
  out <- vMF(sample, method = "Mardia")
  expect_equal(out$m , m, tolerance = 1E-1) #moment estimate part
  cdabyppi:::expect_lt_v(abs(out$k - k), 3 * out$SE$k)
})

test_that("vMF_Full() function works", {
  set.seed(123)
  k <- 3
  m <- c(1, -1) / sqrt(2)
  km <-  k * m
  sample <- movMF::rmovMF(1000, km)
  out <- vMF(sample, method = "smfull")
  cdabyppi:::expect_lt_v(abs(out$km - km), 3 * out$SE$km)

  #with a fixed component
  inkm <- km
  inkm[2] <- NA
  out <- vMF(sample, km = inkm, method = "smfull")
  cdabyppi:::expect_lte_v(abs(out$km - km), 3 * out$SE$km)
})

test_that("vMF() fitting works on dimension 5", {
  set.seed(123)
  p <- 5
  k <- 3
  m <- runif(p, min = -10, 10)
  m <- m / sqrt(sum(m^2))
  km <-  k * m
  sample <- Directional::rvmf(1000, m, k)
  #full method
  out <- vMF(sample, method = "smfull", control = list(tol = 1E-10))
  cdabyppi:::expect_lt_v(abs(out$km - km), 3 * out$SE$km)
  #full with a fixed components
  inkm <- km
  inkm[2] <- NA
  inkm[p] <- NA
  out <- vMF(sample, km = inkm, method = "smfull")
  cdabyppi:::expect_lte_v(abs(out$km - km), 3 * out$SE$km)
  expect_equal(out$SE$km[!is.na(inkm)], rep(0, sum(!is.na(inkm))))

  #Mardia method
  out <- vMF(sample, method = "Mardia")
  expect_equal(out$m , m, tolerance = 1E-1) #moment estimate part
  cdabyppi:::expect_lt_v(abs(out$k - k), 3 * out$SE$k)
})

test_that("vMF matches for simulated weights", {
  set.seed(123)
  p <- 3
  k <- 3
  m <- runif(p, min = -10, 10)
  m <- m / sqrt(sum(m^2))
  km <-  k * m
  Y <- Directional::rvmf(1000, m, k)
  #simulate weights
  set.seed(1342)
  vw <- virtualweights(Y)

  sim1 <- vMF(vw$newY, method = "Mardia")
  dir1 <-  vMF(Y, method = "Mardia", w = vw$w)
  expect_equal(sim1, dir1, tolerance = 1E-3)

  sim2 <- vMF(vw$newY, method = "smfull")
  dir2 <-  vMF(Y, method = "smfull", w = vw$w)
  expect_equal(sim2, dir2, tolerance = 1E-3)
})


test_that("vMF() robust fitting works on dimension 5", {
  set.seed(123)
  p <- 5
  k <- 3
  m <- rep(1, p) #uniform direction poorness due to outliers is evenly distributed in each element
  m <- m / sqrt(sum(m^2))
  km <-  k * m
  set.seed(121)
  sample <- Directional::rvmf(1000, m, k)
  # add outliers
  set.seed(2151)
  outliers <- Directional::rvmf(100, -m, k)
  sample_o <- rbind(sample, outliers)
  #full method, not robust
  out1 <- vMF(sample_o, method = "smfull", control = list(tol = 1E-10))
  #full method, robust, expect to be closer to true value (due to the outliers)
  out2 <- vMF(sample_o, method = "smfull", control = list(tol = 1E-10), cW = 0.1)
  expect_true(all(abs(out2$theta - km) < abs(out1$km - km)))

  #check that fixed components remain fixed for full method
  inkm <- km
  inkm[2] <- NA
  inkm[p] <- NA
  out <- vMF(sample_o, km = inkm, method = "smfull", cW = 0.1)
  expect_equal(out$km[!is.na(inkm)], km[!is.na(inkm)])

  #Mardia method
  out1 <- vMF(sample_o, method = "Mardia")
  #full method, robust, expect to be closer to true value (due to the outliers)
  out2 <- vMF(sample_o, method = "Mardia", cW = 0.1)
  expect_true(all(abs(out2$km - km) < abs(out1$km - km)))
  #expect partially robust Mardia method to not be partially closer
  out3 <- vMF(sample_o, method = "Mardia_robustsm", cW = 0.1)
  expect_equal(out3$m, out1$m) #because m is estimated the same way
  expect_true(all(abs(out3$km - km) > abs(out1$km - km)))
})

test_that("vMF() partially robust fitting works on dimension 5", {
  set.seed(123)
  p <- 5
  k <- 3
  m <- rep(1, p) #uniform direction poorness due to outliers is evenly distributed in each element
  m <- m / sqrt(sum(m^2))
  km <-  k * m
  set.seed(121)
  sample <- Directional::rvmf(1000, m, k)
  # add outliers in concentration only
  set.seed(2151)
  outliers <- Directional::rvmf(100, m, 10 * k)
  sample_o <- rbind(sample, outliers)

  #Mardia method
  out1 <- vMF(sample_o, method = "Mardia")
  #full method, robust, expect to be closer to true value
  out2 <- vMF(sample_o, method = "Mardia", cW = 0.1)
  expect_true(all(abs(out2$km - km) < abs(out1$km - km)))
  #expect partially robust Mardia method to be better at k
  out3 <- vMF(sample_o, method = "Mardia_robustsm", cW = 0.1)
  expect_equal(out3$m, out1$m) #because m is estimated the same way
  expect_lt(abs(out3$k - k), abs(out1$k - k))
})


test_that("controls of FixedPoint() and Rcgmin() are correctly passed", {
  set.seed(123)
  p <- 3
  k <- 3
  m <- runif(p, min = -10, 10)
  m <- m / sqrt(sum(m^2))
  km <-  k * m
  set.seed(123)
  Y <- Directional::rvmf(1000, m, k)

  out_default <- vMF(Y, method = "Mardia", cW = 0.1) #use this packages defaults this is pretty fussy!

  suppressWarnings(out1 <- vMF(Y, method = "Mardia", cW = 0.1,
             control = list(MaxIter = 2, #Fixed point iterations of only 2
                            maxit = 1)))  #Rcgmin iterations of only 1 - warnings of non-convergence
  expect_equal(out1$optim$fpevals, 2)
  expect_error(expect_equal(out_default, out1))

  # expect a different result when Rcgmin package defaults used
  out2 <- vMF(Y, method = "Mardia", cW = 0.1,
           control = list(MaxIter = 2))
  expect_error(expect_equal(out1$km, out2$km))
  expect_error(expect_equal(out_default$km, out2$km))

  # expect a different result when FixedPoint() package defaults used
  suppressWarnings(out3 <- vMF(Y, method = "Mardia", cW = 0.1,
              control = list(maxit = 1)))
  expect_error(expect_equal(out1, out3))
  expect_error(expect_equal(out_default, out3))
})


test_that("dmovMF() and dmvf() are NOT equal", {
  set.seed(123)
  p <- 5
  k <- 3
  m <- runif(p, min = -10, 10)
  m <- m / sqrt(sum(m^2))
  km <-  k * m
  sample <- Directional::rvmf(5, m, k)

  expect_error(expect_equal(movMF::dmovMF(sample, km, log = TRUE),
               Directional::dvmf(sample, k, m, logden = TRUE)))
})

