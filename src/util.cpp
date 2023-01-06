# include "util.hpp"

Rcpp::XPtr< manifold<a1type> > pmanifold(std::string manifoldname){
  manifold<a1type> * out;  //returning a pointer
  if (manifoldname.compare("sphere") == 0){
    out = new mantran::Spos<a1type>();
  } else if (manifoldname.compare("simplex") == 0){
    out = new mantran::simplex<a1type>();
  } else if (manifoldname.compare("Ralr") == 0){
    out = new mantran::Ralr<a1type>();
  } else if (manifoldname.compare("Snative") == 0){
    out = new mantran::Snative<a1type>();
  } else if (manifoldname.compare("Hclr") == 0){
    out = new mantran::Hclr<a1type>();
  } else {
    Rcpp::stop("Manifold not found");
  }

  Rcpp::XPtr< manifold<a1type> > pout(out, true);
  pout.attr("name") = manifoldname;
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

veca1 ptoM(Rcpp::XPtr< manifold<a1type> > pman, veca1 u_ad){
  veca1 z_ad(u_ad.size());
  z_ad = pman->toM(u_ad);
  return(z_ad);
}


vecd pTaylorApprox(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                     vecd u, vecd centre,
                     vecd dynparam, size_t order){
  vecd out(pfun->Range());
  pfun->new_dynamic(dynparam);
  out = taylorapprox(*pfun,
                     centre,
                     order,
                     u);

  return(out);
}


CppAD::ADFun<double> tapefromM(veca1 z,
                               manifold<a1type> *pman){
  //tape relationship between z and h2
  CppAD::Independent(z);
  // range space vector
  veca1 y(0); // vector of ranges space variables - length to be set by from(z);
  y = pman->fromM(z);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(z, y);
  tape.optimize(); //remove some of the extra variables that were used for recording the ADFun, but aren't needed anymore (hopefully very good when logdetJ is constant.
  return(tape);
}


Rcpp::XPtr< CppAD::ADFun<double> > ptapefromM(veca1 z,
                               Rcpp::XPtr<manifold<a1type> > pman){
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  *out = tapefromM(z, pman.checked_get());

  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

