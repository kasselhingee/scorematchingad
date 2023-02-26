# include "exposetransform.h"

// manifold object 'factory'
transform<a1type> * newtransform(const std::string &name){
  transform<a1type> * out;  //returning a pointer
  if (name.compare("alr") == 0){
    out = new mantran::alr<a1type>();
  } else {
    Rcpp::stop("Manifold not found");
  }

  return(out);
}

RCPP_MODULE(transforms) {
  Rcpp::class_< transform_a1type >("transform_ad")
      .factory<const std::string &>(newtransform)
      .method("toM", &transform_a1type::toM, "transform a vector to the manifold")
      .method("fromM", &transform_a1type::fromM, "reverse of toM()")
      .method("logdetJfromM", &transform_a1type::logdetJfromM, "compute the log of the determinant of the Jacobian of fromM()")
  ;
}

