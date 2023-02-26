#ifndef exposemanifold 
#define exposemanifold

#include "scorecompdir_types.h"
#include "mycpp/mantrans.hpp"
#include "mycpp/transforms.hpp"
#include <Rcpp.h>

////////////// Create Pointers to Manifold Objects ///////////////
//in R store a pointer to the ADFun object
//' @noRd
//' @title Generate manifold with transformation object
//' @param manifoldname The name of the manifold to transform to. Either 'sphere' or 'simplex'
//' @return An RCpp::XPtr object pointing to the C++ manifold object
//' @details
//'  + "sphere" for square-root transformation from the simplex to the positive orthant of the sphere
//'  + "simplex" for the simplex without any transformation.
//'  + "Ralr" for the additive log-ratio transformation from the simplex to Euclidean space, using the final component of vectors in the denominator of the ratio.
//'  + "Snative" for the sphere without any transformation
// [[Rcpp::export]]
Rcpp::XPtr< manifold<a1type> > pmanifold(std::string manifoldname);

#endif

