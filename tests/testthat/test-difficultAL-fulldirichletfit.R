test_that("full ppi estimates are with 3 SE for difficult AL with large maxden, p = 5", {
  set.seed(1273)
  p = 5
  ALs <- exp(rsymmetricmatrix(p-1, -4, 4))
  bL <- rep(0, p-1)
  beta <- c(-0.7, -0.8, -0.3, 0, 0)
  set.seed(1111) #this seed leads to some ginormous elements for the second diagonal element of ALs
  suppressMessages(prop <- rhybrid(1000, p, beta, ALs, bL, 35)$samp3) #rhybrid_singly took 1005 seconds, rhybrid() took 13seconds
  #prop %>% as_tibble() %>% tidyr::pivot_longer(everything()) %>% ggplot() + facet_wrap(vars(name)) + geom_freqpoly(aes(x=value))

  est_cppad <- ppi(prop, bL = bL, man = "sphere", bdryweight = "minsq",
                         method = "cppad",
                         acut = 0.01,
                         bdrythreshold = 1E-20,
                         control = list(tol = 1E-10))
  expect_absdiff_lte_v(est_cppad$est$ALs, ALs, 3 * est_cppad$SE$ALs)
  expect_absdiff_lte_v(est_cppad$est$beta, beta, 3 * est_cppad$SE$beta)

  #don't expect that the beta are within a fraction of the true values
  expect_error(expect_equal(est_cppad$est$beta, beta, tolerance = 1E-1))
})

test_that("full ppi estimates are within 3 SE of beta for difficult AL with large maxden, p = 3", {
  set.seed(12735)
  p = 3
  ALs <- exp(rsymmetricmatrix(p-1, -4, 4))
  bL <- rep(0, p-1)
  beta <- c(-0.7, -0.3, 0)
  set.seed(11112) #this seed leads to some ginormous elements for the second diagonal element of ALs
  suppressMessages(prop <- rhybrid(20000, p, beta, ALs, bL, 20)$samp3)
  expect_equal(colMeans(prop), c(0.1, 0.99, 0.1), tolerance = 0.5)

  est_cppad <- ppi(prop, bL = bL, man = "sphere", bdryweight = "minsq",
                         method = "cppad",
                         acut = 0.01,
                         bdrythreshold = 1E-20,
                         control = list(tol = 1E-10))
  expect_error(expect_absdiff_lte_v(est_cppad$est$ALs, ALs, 3 * est_cppad$SE$ALs))
  expect_absdiff_lte_v(est_cppad$est$beta, beta, 3 * est_cppad$SE$beta)
})
