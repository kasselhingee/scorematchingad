#include <scorematchingad_forward.h>
#include <Rcpp.h>

//both Rcpp::wrap and Rcpp::as use copy operations, which ADFun explicitly does not have (the default is 'deleted').
//So exposing pointers of ADFun objects to make interrogation possible
//Hopefully the wrapping for Rcpp::XPtr won't override anything I try here

typedef Rcpp::XPtr<CppAD::ADFun<double>> pADFundouble; //specialisation of ADFun to make exposing to R easier
//RCPP_EXPOSED_WRAP(pADFundouble) //only do wrap which put pointers around class. Rcpp::as() does copying so I will be avoiding it - I'll have to do something bespoke

size_t size_order(pADFundouble * pp){//pointer to a pointer here!
    return (*pp)->size_order();
}

RCPP_MODULE(cppad_module) {
    Rcpp::class_<pADFundouble>("ADFun")
        .property("size_order", &size_order);
}
