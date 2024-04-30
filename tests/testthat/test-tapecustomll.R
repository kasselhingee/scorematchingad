# for interactive testing load_all() with install doesn't include scorematchingad.h properly
test_that("custom ll can be taped", {
  psimplex <- manifoldtransform("sim", "identity", "sim")
  
  lltape <- tapell("
  a1type dirichlet(const veca1 &u, const veca1 &beta) {
          size_t d  = u.size();
          a1type y(0.);  // initialize summation
          for(size_t i = 0; i < d; i++)
          {   y   += beta[i] * log(u[i]);
          }
          return y;
  }", c(0.1, 0.4, 0.5), rep(NA, 3), psimplex$tran, verbose = FALSE)

  expect_equal(pForward0(lltape$ptr, rep(1/3, 3), rep(-0.5, 3)), 3 * (-0.5 * log(1/3)))
  expect_equal(pJacobian(lltape$ptr, rep(1/3, 3), rep(-0.5, 3)), rep(-0.5 * 3, 3))
})

test_that("custom ll errors correctly with wrong signature", {
  psimplex <- manifoldtransform("sim", "identity", "sim")
  expect_error({dirll <- tapell("
a1type dirichlet(const vecd &u, const veca1 &beta) {
        size_t d  = u.size();
        a1type y(0.);  // initialize summation
        for(size_t i = 0; i < d; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
}", c(0.1, 0.4, 0.5), rep(NA, 3), psimplex$tran, verbose = FALSE)},
               "should.*const.*veca1")
})

test_that("score matching with customll matches inbuilt score matching", {
  set.seed(13411)
  mod <- rppi_egmodel(100)
  tran <- manifoldtransform("sim", "sqrt", "sph")$tran
  lltape <- tapell("
  a1type customppi(const veca1 &u, const veca1 &theta) {
          a1type y;
          y = ll::ll_ppi(u, theta);
          return y;
  }", c(0.1, 0.4, 0.5), ppi_paramvec(p = 3, betap = tail(mod$beta, 1)), tran)

  Y <- mod$sample

  tapes_custom <- buildsmdtape(
     start = "sim",
     tran = "alr",
     end = "Euc",
     ll = "
  a1type customppi(const veca1 &u, const veca1 &theta) {
          a1type y;
          y = ll::ll_ppi(u, theta);
          return y;
  }",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3, betap = tail(mod$beta, 1)),
     bdryw = "ones",
     verbose = FALSE)
  tapes <- buildsmdtape(
     start = "sim",
     tran = "alr",
     end = "Euc",
     ll = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3, betap = tail(mod$beta, 1)),
     bdryw = "ones",
     verbose = FALSE)

  est_custom <- cppad_closed(tapes_custom$smdtape, Y) 
  est <- cppad_closed(tapes$smdtape, Y) 
  expect_equal(est_custom$est, est$est)
  expect_equal(est_custom$covar, est$covar)
})
