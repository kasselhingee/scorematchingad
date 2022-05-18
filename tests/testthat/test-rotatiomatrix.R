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
