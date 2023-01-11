#include "scorecompdir_types.h"
#include "mycpp/mantrans.hpp"
#include <Rcpp.h>
using namespace Rcpp;

// manifold object 'factory'
manifold<a1type> * newmanifold(const std::string &manifoldname){
  manifold<a1type> * out;  //returning a pointer
  if (manifoldname.compare("sphere") == 0){
    out = new mantran::Spos<a1type>();
  } else if (manifoldname.compare("simplex") == 0){
    out = new mantran::simplex<a1type>();
  } else if (manifoldname.compare("Ralr") == 0){
    out = new mantran::Ralr<a1type>();
  } else if (manifoldname.compare("Snative") == 0){
    out = new mantran::Snative<a1type>();
  } else {
    Rcpp::stop("Manifold not found");
  }

  return(out);
}

RCPP_MODULE(manif) {
  Rcpp::class_< manifold<a1type> >("manifold")
      .factory<const std::string &>(newmanifold)
      .method("toM", &manifold<a1type>::toM);

}
