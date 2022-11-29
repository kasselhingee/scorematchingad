skip_on_cran()
test_that("Sphere manifold passes lightweight standard tests",{
  pman <- manifoldtransform("sphere")
  u <- c(0.1, 0.2, 1 - 0.3)
  out <- testmanifold(pman, u)
  expect_equal(out, 0)
})

test_that("Simplex manifold passes lightweight standard tests",{
  pman <- manifoldtransform("simplex")
  u <- c(0.1, 0.2, 1 - 0.3)
  out <- testmanifold(pman, u)
  expect_equal(out, 0)
})

test_that("Ralr manifold passes lightweight standard tests",{
  pman <- manifoldtransform("Ralr")
  u <- c(0.1, 0.2, 1 - 0.3)
  out <- testmanifold(pman, u)
  expect_equal(out, 0)
})

test_that("Snative manifold object passes lightweight standard tests", {
  pman <- manifoldtransform("Snative")
  u <- c(0.1, 0.2)
  u <- u / sqrt(sum(u^2))
  out <- testmanifold(pman, u)
  expect_equal(out, 0)
})

