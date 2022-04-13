#include <Rcpp.h>
using namespace Rcpp ;
# include <cppad/cppad.hpp>
typedef CppAD::ADFun<double> ADFund;

RCPP_EXPOSED_CLASS_NODECL(class_ADFun)
RCPP_MODULE(class_ADFun) {
    class_<ADFund>("ADFun")
    .constructor()
    .method("Jacobian", &ADFund::Jacobian)
    ;
}
