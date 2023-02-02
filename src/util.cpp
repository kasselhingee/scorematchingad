# include "util.h"

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

