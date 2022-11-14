# ifndef ADFUN_MODULE
# define ADFUN_MODULE

// Exposing some of the CppAD ADFun class, with value of type double
// then RCPP_EXPOSED_AS and RCPP_EXPOSED_WRAP can be used to export or read in the objects
// see cppad/core/ad_fun.hpp for the definition of ADFun in CppAD


#include "scorecompdir_types.h"
//#include "scm_interface.cpp"
# include <cppad/cppad.hpp>
# include <Rcpp.h>
using namespace Rcpp;

// class ADFun_double defined as a special version of the templated ADFun
// following how RcppAnnoy has done it
RCPP_EXPOSED_CLASS_NODECL(ADFun_double)
RCPP_MODULE(ADFun_double) {

  class_<ADFun_double>("ADFun")

  .constructor() //the default constructor

  //would love the copy constructor too, but want to stop R from actually copying

  //.method("Domain", &CppAD::ADFun<double>::Domain, "The dimension of the domain (i.e. the length of the vector that this ADFun evaluates")

  ; //in the Rcpp vignette on modules there is a single semi-colon at the end of the definitions
}


# endif
