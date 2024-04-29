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
