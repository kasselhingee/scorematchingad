#include "scorecompdir_types.h"
#include "mycpp/mantrans.hpp"
#include <Rcpp.h>
using namespace Rcpp;

// manifold object 'factory'
manifold<a1type> * newmantran(const std::string &manifoldname){
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

manifold<a1type> * newmantran(const manifold<a1type> &inmantran){
  manifold<a1type> * out = newmantran(inmantran.name());
  return(out);
}

typedef manifold<a1type> manifold_a1type;

RCPP_MODULE(manifolds) {
  Rcpp::class_< manifold_a1type >("mantran_ad")
      .factory<const std::string &>(newmantran)
      //.factory<const manifold_a1type &>(newmantran)
      .method("toM", &manifold_a1type::toM, "transform a vector to the manifold")
      .method("fromM", &manifold_a1type::fromM, "reverse of toM()")
      .method("logdetJfromM", &manifold_a1type::logdetJfromM, "compute the log of the determinant of the Jacobian of fromM()")
      .method("Pmatfun", &manifold_a1type::Pmatfun, "Pmatfun(z) returns the matrix that orthogonally projects onto the manifold's tangent space at z")
      .method("dPmatfun", &manifold_a1type::dPmatfun, "dPmatfun(z, i) returns the element-wise derivative of Pmatfun() at location z with respect to the ith dimension")
  ;
}

//adadsafRCPP_EXPOSED_CLASS_NODECL(manifold_a1type) //so that Rcpp knows how to wrap and unwrap
