skip_on_cran() #not relevant anymore - established issue was simulation via Directional::rvmf() simulation
test_that("vMF objective, grad, hess calculations are repeatable", {
  set.seed(123)
  p <- 3
  k <- 3
  m <- runif(p, min = -10, 10)
  m <- m / sqrt(sum(m^2))
  theta <- k*m
  set.seed(12432)
  sample <- movMF::rmovMF(1000, theta)
  thetatape <- c(1, -1, 1)

  p <- length(theta)
  tapes <- buildsmotape("Snative", "vMF",
                        rep(1, p)/sqrt(p), rep(NA, p),
                        weightname = "ones",
                        verbose = FALSE)
  # for true theta
  obj <- replicate(1000, smobj(tapes$smotape, theta, sample))
  expect_lt(sum(range(obj) * c(1, -1)), 1E-10)
  grads <- replicate(1000, smobjgrad(tapes$smotape, theta, sample))
  cdabyppi:::expect_lt_v(apply(grads, MARGIN = 1, function(x){sum(range(x) * c(1, -1))}), rep(1E-10, nrow(grads)))
  hesss <- replicate(1000, smobjhess(tapes$smotape, theta, sample))
  cdabyppi:::expect_lt_v(apply(hesss, MARGIN = c(1,2), function(x){sum(range(x) * c(1, -1))}), rep(1E-10, p^2))

  # for incorrect theta
  obj <- replicate(1000, smobj(tapes$smotape, 2*thetatape, sample))
  expect_lt(sum(range(obj) * c(1, -1)), 1E-10)
  grads <- replicate(1000, smobjgrad(tapes$smotape, 2*thetatape, sample))
  cdabyppi:::expect_lt_v(apply(grads, MARGIN = 1, function(x){sum(range(x) * c(1, -1))}), rep(1E-10, nrow(grads)))
  hesss <- replicate(1000, smobjhess(tapes$smotape, 2*thetatape, sample))
  cdabyppi:::expect_lt_v(apply(hesss, MARGIN = c(1,2), function(x){sum(range(x) * c(1, -1))}), rep(1E-10, p^2))
})

test_that("vMF cppadest is repeatable for fixed tape", {
  set.seed(123)
  p <- 3
  k <- 3
  m <- runif(p, min = -10, 10)
  m <- m / sqrt(sum(m^2))
  theta <- k*m
  set.seed(12432)
  sample <- movMF::rmovMF(1000, theta)
  thetatape <- c(1, -1, 1)

  p <- length(theta)
  tapes <- buildsmotape("Snative", "vMF",
                        rep(1, p)/sqrt(p), rep(NA, p),
                        weightname = "ones",
                        verbose = FALSE)
  kests <- replicate(100, sqrt(sum(cppadest(tapes$smotape, thetatape, sample)$par^2)))
  expect_lt(sum(range(kests) * c(1, -1)), 1E-10)
})

test_that("vMF different tape objects have same estimates", {
  set.seed(123)
  p <- 3
  k <- 3
  m <- runif(p, min = -10, 10)
  m <- m / sqrt(sum(m^2))
  theta <- k*m
  set.seed(12432)
  sample <- movMF::rmovMF(1000, theta)
  thetatape <- c(1, -1, 1)
  p <- length(theta)

  tapenest <- function(){
    tapes <- buildsmotape("Snative", "vMF",
                          rep(1, p)/sqrt(p), rep(NA, p),
                          weightname = "ones",
                          verbose = FALSE)
    sqrt(sum(cppadest(tapes$smotape, thetatape, sample)$par^2))
  }
  kests <- replicate(100, tapenest())
  expect_lt(sum(range(kests) * c(1, -1)), 1E-10)
})

test_that("Estimating repeated vMF estimation for the same data gives identical results, even when Rcgmin iterations are low", {
  set.seed(123)
  p <- 3
  k <- 3
  m <- runif(p, min = -10, 10)
  m <- m / sqrt(sum(m^2))
  km <-  k * m
  set.seed(1241)
  Y <- movMF::rmovMF(100, m * k)

  kests <- replicate(100, {
    sim1 <- vMF(Y, method = "Mardia")
    sim1$est$k })
  expect_lt(sum(range(kests) * c(1, -1)), 1E-10)
})
