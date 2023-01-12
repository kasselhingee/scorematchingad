# include "exposemanifold.hpp"

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

Rcpp::XPtr< manifold<a1type> > pmanifold(std::string manifoldname){
  manifold<a1type> * out;  //returning a pointer
  out = newmantran(manifoldname);
  Rcpp::XPtr< manifold<a1type> > pout(out, true);
  pout.attr("name") = out->name();
  return(pout);
}

int testmanifold(Rcpp::XPtr< manifold<a1type> > pman, veca1 u_ad){
  Rcpp::Rcout << "Starting tests" << std::endl;
  // toM then fromM get back to u
  Rcpp::Rcout << "               Input u was: " << u_ad.transpose() << std::endl;
  veca1 z_ad(u_ad.size());
  z_ad = pman->toM(u_ad);
  Rcpp::Rcout << "                 After toM: " << z_ad.transpose() << std::endl;
  veca1 u2_ad(u_ad.size());
  u2_ad = pman->fromM(z_ad);
  Rcpp::Rcout << "      After toM then fromM: " << u2_ad.transpose() << std::endl;
  if ((u2_ad - u_ad).array().abs().maxCoeff() > 1E-8){
    Rcpp::Rcout << "toM then fromM not passed." << std::endl;
    return(1);
  }

  // Run the other elements
  Rcpp::Rcout << " logdetJ_fromM at toM(u): " << pman->logdetJfromM(z_ad) << std::endl;
  Rcpp::Rcout << " Pmat at toM(u): " << std::endl << pman->Pmatfun(z_ad) << std::endl;
  for (long int d=0; d<u_ad.size(); d++){
    Rcpp::Rcout << " dPmat at toM(u) in dimension " << d <<":" << std::endl << pman->dPmatfun(z_ad, d) << std::endl;
  }
  return(0);
}

RCPP_MODULE(manifolds) {
  Rcpp::class_< manifold_a1type >("mantran_ad")
      .factory<const std::string &>(newmantran)
      .method("toM", &manifold_a1type::toM, "transform a vector to the manifold")
      .method("fromM", &manifold_a1type::fromM, "reverse of toM()")
      .method("logdetJfromM", &manifold_a1type::logdetJfromM, "compute the log of the determinant of the Jacobian of fromM()")
      .method("Pmatfun", &manifold_a1type::Pmatfun, "Pmatfun(z) returns the matrix that orthogonally projects onto the manifold's tangent space at z")
      .method("dPmatfun", &manifold_a1type::dPmatfun, "dPmatfun(z, i) returns the element-wise derivative of Pmatfun() at location z with respect to the ith dimension")
  ;
}

