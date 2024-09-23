#include <scorematchingad_forward.h>
#include <Rcpp.h>

RCPP_EXPOSED_WRAP(ADFundouble) //only do wrap which put pointers around class. Rcpp::as() does copying so I will be avoiding it - I'll have to do something bespoke

RCPP_MODULE(cppad_module) {
    using namespace Rcpp;

    class_<ADFundouble>("ADFun")
        .method("size_order", &ADFundouble::size_order);
}
