#include "scorematchingad_forward.h"
#include <Rcpp.h>

// Expose the ADFun class and its methods
RCPP_MODULE(cppad_module) {
    using namespace Rcpp;

    class_<ADFundouble>("ADFun")
        .constructor<>()
        //.method("Forward", (std::vector<double> (ADFundouble::*)(size_t, const std::vector<double>&, std::ostream&)) &ADFundouble::Forward)
        //.method("Reverse", &ADFundouble::Reverse)
        //.method("Jacobian", &ADFundouble::Jacobian)
        //.method("Hessian", &ADFundouble::Hessian)
        //.method("optimize", &ADFundouble::optimize)
        //.method("size_var", &ADFundouble::size_var)
        .method("size_order", &ADFundouble::size_order);
}

