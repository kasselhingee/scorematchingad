test_that("getting inline plugin passes", {
  expect_no_error(inline::getPlugin("scorematchingad"))
})

test_that("simple cxx compiles and runs", {
  expect_no_error({fx <- inline::cxxfunction(signature(x = "integer", y = "numeric"), 
                             "return wrap( as<int>(x) * as<double>(y));",
                             plugin = "scorematchingad" )})
  expect_no_error(fx(2L, 5))
  #write(strsplit(fx@code, "\n")[[1]], file = "tmp.cpp")
})

test_that("Rcpp::cppFunction() accesses Eigen and wraps automatically", {
# Rcpp::cppFunction() allows for automatic use of wrap and as!
# note it is using the inline plugin
fxptr <- Rcpp::cppFunction("
double foo(vecd x) {
 double out;
 out = x.sum();
 return out;
}
", depends = c("RcppEigen", "scorematchingad"))
  expect_equal(fxptr(c(1, 2)), 3)
})


test_that("RcppXPtrUtils accesses Eigen and wraps automatically", {
# Rcpp::cppFunction() allows for automatic use of wrap and as!
# note it is using the inline plugin
{ptr <- RcppXPtrUtils::cppXPtr("
double foo(vecd x) {
 double out;
 out = x.sum();
 return out;
}
", depends = c("RcppEigen", "scorematchingad"), verbose = TRUE, showOutput = FALSE)} |> expect_no_error()
  # can't checkXPtr because argument type has commas in it
})


Rcpp::cppFunction("
veca1 foo(transform<a1type> & tran, veca1 z) {
 veca1 u(0);
 //u = tran.fromM(z);
 return u;
}
", depends = c("RcppEigen", "scorematchingad"))


#RcppXPtrUtils looks like exactly what I want

test_that("returning pointer can be used elsewhere", {
  expect_no_error({
fx <- inline::cxxfunction(
signature(z_ad = "numeric", theta_ad = "numeric", llname = "string", 
body = "

a1type (*ll)(const veca1 &, const veca1 &) = nullptr;
ll = ll::ll_ppi;
return Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
",
                             plugin = "scorematchingad" )
})
  expect_no_error(fx(2L, 5))
  #write(strsplit(fx@code, "\n")[[1]], file = "tmp.cpp")
})
