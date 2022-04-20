test_that("On sphere with prodsq, objective values independent of tape", {
  u = matrix(c(0.001, 0.011, 1 - 0.01 - 0.011), nrow = 1)

  tapeval1 = c(u,3,3,3)
  smofun1 <- ptapesmo(tapeval1, 3,
                     manifoldname = "sphere", "prodsq",
                     acut = 0.1)

  tapeval2 = c(0.2, 0.2, 0.6,3,3,3)
  smofun2 <- ptapesmo(tapeval2, 3,
                      manifoldname = "sphere", "prodsq",
                      acut = 0.1)

  tapeval3 = c(u[c(2, 1, 3)],3,3,3)
  smofun3 <- ptapesmo(tapeval3, 3,
                      manifoldname = "sphere", "prodsq",
                      acut = 0.1)

  expect_equal(smobj(smofun1, c(1, 1, 1), u), smobj(smofun2, c(1, 1, 1), u))
  expect_equal(smobj(smofun1, c(1, 1, 1), u), smobj(smofun3, c(1, 1, 1), u))
  expect_equal(smobjgrad(smofun1, c(1, 1, 1), u), smobjgrad(smofun2, c(1, 1, 1), u))
  expect_equal(smobjgrad(smofun1, c(1, 1, 1), u), smobjgrad(smofun3, c(1, 1, 1), u))
})

test_that("On sphere with minsq, objective values independent of tape", {
  u = matrix(c(0.001, 0.011, 1 - 0.01 - 0.011), nrow = 1)

  tapeval1 = c(u,3,3,3)
  smofun1 <- ptapesmo(tapeval1, 3,
                      manifoldname = "sphere", "minsq",
                      acut = 0.1)

  tapeval2 = c(0.2, 0.2, 0.6,3,3,3)
  smofun2 <- ptapesmo(tapeval2, 3,
                      manifoldname = "sphere", "minsq",
                      acut = 0.1)

  tapeval3 = c(u[c(2, 1, 3)],3,3,3)
  smofun3 <- ptapesmo(tapeval3, 3,
                      manifoldname = "sphere", "minsq",
                      acut = 0.1)

  expect_equal(smobj(smofun1, c(1, 1, 1), u), smobj(smofun2, c(1, 1, 1), u))
  expect_equal(smobj(smofun1, c(1, 1, 1), u), smobj(smofun3, c(1, 1, 1), u))
  expect_equal(smobjgrad(smofun1, c(1, 1, 1), u), smobjgrad(smofun2, c(1, 1, 1), u))
  expect_equal(smobjgrad(smofun1, c(1, 1, 1), u), smobjgrad(smofun3, c(1, 1, 1), u))
})
