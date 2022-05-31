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
