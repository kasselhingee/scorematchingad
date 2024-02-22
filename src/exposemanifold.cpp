# include "exposemanifold.h"

// manifold object 'factory'
manifold<a1type> * newmantran(const std::string &manifoldname){
  manifold<a1type> * out;  //returning a pointer
  if (manifoldname.compare("sph") == 0){
    out = new mantran::sph<a1type>();
  } else if (manifoldname.compare("sim") == 0){
    out = new mantran::sim<a1type>();
  } else if (manifoldname.compare("Euc") == 0){
    out = new mantran::Euc<a1type>();
  } else if (manifoldname.compare("Hn111") == 0){
    out = new mantran::Hn111<a1type>();
  } else {
    Rcpp::stop("Manifold not found");
  }

  return(out);
}

// transform factory
transform<a1type> * newtransform(const std::string &name){
  transform<a1type> * out;  //returning a pointer
  if (name.compare("alr") == 0){
    out = new mantran::alr<a1type>();
  } else if (name.compare("clr") == 0){
    out = new mantran::clr<a1type>();
  } else if (name.compare("sqrt") == 0){
    out = new mantran::sqrt<a1type>();
  } else if (name.compare("identity") == 0){
    out = new mantran::identity<a1type>();
  } else {
    Rcpp::stop("Transform not found");
  }

  return(out);
}

RCPP_MODULE(manifolds) {
  Rcpp::class_< manifold_a1type >("man_ad")
      .factory<const std::string &>(newmantran)
      .method("Pmatfun", &manifold_a1type::Pmatfun, "Pmatfun(z) returns the matrix that orthogonally projects onto the manifold's tangent space at z")
      .method("dPmatfun", &manifold_a1type::dPmatfun, "dPmatfun(z, i) returns the element-wise derivative of Pmatfun() at location z with respect to the ith dimension")
      .method("name", &manifold_a1type::name)
  ;
  
  Rcpp::class_< transform_a1type >("transform_ad")
      .factory<const std::string &>(newtransform)
      .method("toM", &transform_a1type::toM, "transform a vector to the manifold")
      .method("fromM", &transform_a1type::fromM, "reverse of toM()")
      .method("logdetJfromM", &transform_a1type::logdetJfromM, "compute the log of the determinant of the Jacobian of fromM()")
      .method("name", &transform_a1type::name)
  ;
}

