test_that("Kent distribution is fitted correctly", {
  #simulate
  set.seed(3121)
  p = 3
  Amat_raw <- rsymmetricmatrix(p, -10, 10)
  # zero out an eigen value, the first one
  A_es <- eigen(Amat_raw)
  Amat <- Amat_raw - A_es$values[[1]] * diag(1, nrow(Amat_raw), ncol(Amat_raw))
  mu <- A_es$vectors[, 1]
  A_es <- eigen(Amat)
  stopifnot(any(abs(A_es$values) < 1E-5))
  k = 2
  sample <- rFB(100, k, mu, Amat)


})
