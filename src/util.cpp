# include "util.hpp"

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
