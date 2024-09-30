test_that("fixdynamic() reduces the dimension of the dynamic parameters", {
  ppitape <- tape_uld_inbuilt(name = "ppi",
                  x = c(0.2, 0.3, 0.5),
                  theta = ppi_paramvec(p = 3, AL = 1.1, bL = 1.2, beta = -0.1)) 
  fixedtape <- fixdynamic(ppitape, ppitape$dyntape, t_u2i(ppi_paramvec(p = 3, AL = TRUE, beta = TRUE)))
  expect_equal(fixedtape$size_dyn_ind, 2)
})

