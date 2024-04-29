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
  expect_equal(evalll(dirll, rep(1/3, 3), rep(-0.5, 3)), 3 * (-0.5 * log(1/3)))
  
  psimplex <- manifoldtransform("sim", "identity", "sim")
  lltape <- tapell(dirll, c(0.1, 0.4, 0.5), rep(NA, 3), psimplex$tran, verbose = FALSE)

  expect_equal(pForward0(lltape$ptr, rep(1/3, 3), rep(-0.5, 3)), 3 * (-0.5 * log(1/3)))
  expect_equal(pJacobian(lltape$ptr, rep(1/3, 3), rep(-0.5, 3)), rep(-0.5 * 3, 3))
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

test_that("cppFunction that does tapellcpp() of a custom ll works on Windows", {
  myfun <- Rcpp::cppFunction(
depends =  c("RcppEigen", "scorematchingad", "Rcpp"),
includes = c("#include <utils/PrintFor.hpp>", "#include <tapellcpp.h>", "#include <manifoldtransforms/transforms.hpp>", "#include <utils/wrapas.hpp>"),
verbose = TRUE,
code = "
  Rcpp::XPtr< CppAD::ADFun<double> > tapedirichlet(veca1 z, veca1 theta, transform_a1type & tran, Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta, bool verbose);

  a1type dirichlet(const veca1 &u, const veca1 &beta) {
          PrintForVec(\"u is:\", u);
          size_t d  = u.size();
          a1type y(0.);  // initialize summation
          for(size_t i = 0; i < d; i++)
          {   y   += beta[i] * log(u[i]);
          }
          return y;
  }

  Rcpp::XPtr< CppAD::ADFun<double> > tapedirichlet(veca1 z, veca1 theta, transform_a1type & tran, Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta, bool verbose){
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>;
  *out = tapellcpp(z,
                theta,
                *dirichlet,
                tran,
                fixedtheta,
                verbose);

  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
  }
  ")

  maninfo <- manifoldtransform("sim", "identity", "sim")
  mytape <- myfun(c(0.1, 0.4, 0.5), c(1,1,1), maninfo$tran, rep(FALSE, 3), verbose = TRUE)

  expect_equal(pForward0(mytape, rep(1/3, 3), rep(-0.5, 3)), 3 * (-0.5 * log(1/3)))
  expect_equal(pJacobian(mytape, rep(1/3, 3), rep(-0.5, 3)), rep(-0.5 * 3, 3))

})
