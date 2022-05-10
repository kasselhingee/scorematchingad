
test_that("Sphere manifold passes lightweight standard tests",{
  pman <- pmanifold("sphere")
  u <- c(0.1, 0.2, 1 - 0.3)
  out <- testmanifold(pman, u)
  expect_equal(out, 0)
})

test_that("Ralr manifold passes lightweight standard tests",{
  pman <- pmanifold("Ralr")
  u <- c(0.1, 0.2, 1 - 0.3)
  out <- testmanifold(pman, u)
  expect_equal(out, 0)
})

