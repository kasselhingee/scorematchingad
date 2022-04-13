Rcpp::exposeClass("ADFun",
                  CppClass = "CppAD::ADFun<double>",
                  constructors = NULL, #hopefully just the default (empty) constructor
                  fields = NULL,
                  methods = "Jacobian",
                  header = "# include <cppad/cppad.hpp>")
