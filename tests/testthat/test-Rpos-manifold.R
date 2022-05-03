
test_that("Rpos manifold passes standard tests",{
  Rpos <- pmanifold("Rpos")
  u <- c(0.1, 0.2, 1 - 0.3)
  out <- testmanifold(Rpos, u)
  expect_equal(out, 0)
})

