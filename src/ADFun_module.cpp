#include <scorematchingad_forward.h>
#include <Rcpp.h>


RCPP_MODULE(cppad_module) {
    using namespace Rcpp;

    class_<ADFundouble>("ADFun")
        .method("size_order", &ADFundouble::size_order);
}
