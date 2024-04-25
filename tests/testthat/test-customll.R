# for interactive testing load_all() with install doesn't include scorematchingad.h properly
test_that("customll can compile a ll function that is evaluated by evalll and taping later", {
  dirll <- customll("
  a1type dirichlet(const veca1 &u, const veca1 &beta) {
          size_t d  = u.size();
          a1type y(0.);  // initialize summation
          for(size_t i = 0; i < d; i++)
          {   y   += beta[i] * log(u[i]);
          }
          return y;
  }")
  expect_equal(3 * (-0.5 * log(1/3)), evalll(dirll, rep(1/3, 3), rep(-0.5, 3)))
  
  psimplex <- manifoldtransform("sim", "identity", "sim")
  lltape <- tapell(dirll, rep(1/3, 3), rep(NA, 3), psimplex$tran, verbose = FALSE)

  expect_equal(3 * (-0.5 * log(1/3)), pForward0(lltape, rep(1/3, 3), rep(-0.5, 3)))
  expect_equal(rep(-0.5 * 3, 3), pJacobian(lltape, rep(1/3, 3), rep(-0.5, 3)))
})

test_that("customll errors correctly with wrong signature", {
  expect_error({dirll <- customll("
a1type dirichlet(const vecd &u, const veca1 &beta) {
        size_t d  = u.size();
        a1type y(0.);  // initialize summation
        for(size_t i = 0; i < d; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
}")})
})

