# include "tapell.h"

Rcpp::XPtr< CppAD::ADFun<double> > ptapell2(veca1 z_ad, //data measurement on the M manifold
                                     veca1 theta_ad,
                                     Rcpp::XPtr<llPtr> llfXPtr, //the log likelihood function
                                     transform_a1type & tran,
                                     Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta,
                                     bool verbose
                                     ){
  // unwrap likelihood function
  llPtr func = *Rcpp::XPtr<llPtr>(llfXPtr);

  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  *out = tapellcpp(z_ad,
                theta_ad,
                func,
                tran,
                fixedtheta,
                verbose);

  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}


Rcpp::XPtr<llPtr> getllptr(std::string llname){
  //choose ll function
  llPtr ll = nullptr;
  if (llname.compare("dirichlet") == 0){
    ll = ll::ll_dirichlet;
  }
  if (llname.compare("ppi") == 0){
    ll = ll::ll_ppi;
  }
  if (llname.compare("vMF") == 0){
    ll = ll::ll_vMF;
  }
  if (llname.compare("Bingham") == 0){
    ll = ll::ll_Bingham;
  }
  if (llname.compare("FB") == 0){
    ll = ll::ll_FB;
  }
  //check ll function
  if (ll == nullptr){
    throw std::invalid_argument("Matching ll function not found");
  }

  Rcpp::XPtr< llPtr > pout(new llPtr(ll), true);
  return(pout);
}

 
a1type evalll(Rcpp::XPtr<llPtr> llfXPtr, const veca1& u, const veca1& theta){
  llPtr func = *Rcpp::XPtr<llPtr>(llfXPtr);
  a1type out;
  out = func(u, theta); //implicit dereferencing of function pointer as per: https://www.learncpp.com/cpp-tutorial/function-pointers/
  return(out);
}

