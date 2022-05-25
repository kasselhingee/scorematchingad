test_that("Rotation matrix generation on high dimensions", {
  u <- c(1, 2, 3, 4, 5, 6)
  Rmat <- vec2northpole(u)
  newu <- Rmat %*% u
  expect_equal(newu[1], sqrt(sum(u^2)))
  expect_equal(newu[-1], rep(0, length(u) - 1))

  u2 <- c(10, 100, 3, -10, 7, 0)
  expect_equal(sqrt(sum((Rmat %*% u2)^2)), sqrt(sum(u2^2)))

  u3 <- c(2, 0, 0, 0, 0, 0)
  expect_equal(sqrt(sum((Rmat %*% u3)^2)), 2)

  u4 <- c(6, 5, 4, 3, 2, 1)
  expect_equal(sqrt(sum((Rmat %*% u4)^2)), sqrt(sum(u4^2)))
})

test_that("Rotation matrix with inverse gets to any direction", {
  u <- c(1, 2, 3, 4, 5, 6)
  uend <- c(3, 4, 1, 2, 5, 6)
  Rmat1 <- vec2northpole(u)
  Rmat2 <- vec2northpole(uend)
  Rmat <- solve(Rmat2) %*% Rmat1

  expect_equal(uend, Rmat %*% u, ignore_attr = TRUE) #because size of u == size of uend
})

test_that("Rotation matrix matches Directional::rotation", {
  u <- c(1, 2, 3, 4, 5, 6)
  nthpole <- c(1, 0, 0 , 0, 0 , 0)

  # to northpole
  Rmat1 <- vec2northpole(u/sqrt(sum(u^2)))
  Rmat1B <- Directional::rotation(u/sqrt(sum(u^2)), nthpole)
  expect_equal(Rmat1 %*% u, Rmat1B %*% u, tolerance = 1E-5)
  expect_equal(Rmat1B, Rmat1, tolerance = 1E-5)
  round(Rmat1 %*% solve(Rmat1B), 2)
  round(Rmat1, 2)
  round(Rmat1B, 2)
  round(Rmat1 %*% t(Rmat1), 2)
  expect_equal(Rmat1 %*% solve(Rmat1B), diag(1, nrow = length(u)), tolerance = 1E-1)
  round(Rmat1 %*% solve(Rmat1B) - diag(1, nrow = length(u)), 3)

  uend <- c(3, 4, 1, 2, 5, 6)
  Rmat2 <- vec2northpole(uend / sqrt(sum(uend^2)))
  Rmat <- solve(Rmat2) %*% Rmat1
  RmatB <- Directional::rotation(u/sqrt(sum(u^2)), uend / sqrt(sum(uend^2)))
  round(Rmat, 2)
  round(RmatB, 2)
  expect_equal(Rmat, RmatB, tolerance = 0.5)
  expect_equal(uend, Rmat %*% u, ignore_attr = TRUE) #because size of u == size of uend
  expect_equal(RmatB %*% u, Rmat %*% u)

  utest <- runif(6, -1, 1)
  expect_equal(sqrt(sum((Rmat %*% utest)^2)), sqrt(sum(utest^2)))
  expect_equal(sqrt(sum((RmatB %*% utest)^2)), sqrt(sum(utest^2)))
  expect_equal(RmatB %*% utest, Rmat %*% utest)
})
