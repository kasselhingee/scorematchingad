test_that("fixdynamic() reduces the dimension of the dynamic parameters", {
  ppitape <- tape_uld_inbuilt(name = "ppi",
                  x = c(0.2, 0.3, 0.5),
                  theta = ppi_paramvec(p = 3, AL = 1.1, bL = 1.2, beta = -0.1)) 
  fixedtape <- fixdynamic(ppitape, ppitape$dyntape, t_u2i(ppi_paramvec(p = 3, AL = TRUE, beta = TRUE)))
  expect_equal(fixedtape$size_dyn_ind, 2)
  expect_equal(fixedtape$eval(ppitape$xtape, c(1.2, 1.2)), ppitape$eval(ppitape$xtape, ppitape$dyntape))
  expect_equal(fixedtape$eval(ppitape$xtape, c(0, 0)), ppitape$eval(ppitape$xtape, ppi_paramvec(p = 3, AL = 1.1, bL = 0, beta = -0.1)))
})

test_that("fixindependent() reduces the dimension of the independent arguments", {
  ppitape <- tape_uld_inbuilt(name = "ppi",
                              x = c(0.2, 0.3, 0.5),
                              theta = ppi_paramvec(p = 3, AL = 1.1, bL = 1.2, beta = -0.1)) 
  fixedtape <- fixindependent(ppitape, ppitape$xtape, c(TRUE, FALSE, FALSE))
  expect_equal(fixedtape$domain, 2)
  expect_equal(fixedtape$xtape, c(0.3, 0.5))
  expect_equal(fixedtape$forward(0, c(0.3, 0.5)), ppitape$forward(0, c(0.2, 0.3, 0.5)))
  expect_equal(fixedtape$forward(0, c(0.4, 0.4)), ppitape$forward(0, c(0.2, 0.4, 0.4)))
})

test_that("keeprange() can drop range elements", {
  ppitape <- tape_uld_inbuilt(name = "ppi",
                              x = c(0.2, 0.3, 0.5),
                              theta = ppi_paramvec(p = 3, AL = 1.1, bL = 1.2, beta = -0.1))
  jactape <- tape_Jacobian(ppitape)
  dtape <- keeprange(jactape, c(2, 3))
  expect_equal(dtape$forward(0, dtape$xtape), jactape$forward(0, jactape$xtape)[-1])

  # error when try to remove all the range
  expect_error(keeprange(jactape, vector("integer", 0)))
})
