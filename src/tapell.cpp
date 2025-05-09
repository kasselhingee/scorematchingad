# include "tapell.h"

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
    Rcpp::stop("Matching ll function not found");
  }

  Rcpp::XPtr< llPtr > pout(new llPtr(ll), true);
  return(pout);
}

pADFun tape_uld_inbuilt(std::string name, veca1 x, veca1 theta){
  Rcpp::XPtr < llPtr > ptr = getllptr(name); 
  llPtr func = *ptr;
  CppAD::ADFun<double> tape;
  CppAD::Independent(x, theta);
  veca1 y(1);
  y(0) = func(x, theta);
  tape.Dependent(x, y);
  tape.check_for_nan(false);
  pADFun out(tape, x, theta, name);
  return(out);
}

 
