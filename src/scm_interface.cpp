//for content that is Rcpp specific

#include "mycpp/scm.cpp"
using namespace Rcpp;

//in R store a pointer to the ADFun object
//' @title The score matching objective calculator.
//' @param xbetain a concatenated vector of sqrt(x) and beta
//' @param n The dimension of x.
//' @param manifoldname The name of the manifold to transform to
//' @param weightname The name of the weight function to use
//' @param acut The constraint a_c in the weight function
//' @return An RCpp::XPtr object pointing to the ADFun
//' @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> > ptapesmo(svecd xbetain,
                                      size_t n,
                                      std::string llname,
                                      std::string manifoldname,
                                      std::string weightname,
                                      const double acut){
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer

  //instantiate the Spos manifold
  manifold<a1type> Spos = {
    Spos::toS, Spos::Pmat_S, Spos::dPmat_S,
    Spos::fromS, Spos::logdetJ_fromS,
    };

  manifold<a1type> simplex = {
    simplex::toM, simplex::Pmat_M, simplex::dPmat_M,
    simplex::fromM, simplex::logdetJ_fromM
  };

  //choose weight function
  a1type (*h2fun)(const veca1 &, const double &) = nullptr;
  if (weightname.compare("prodsq") == 0){
    h2fun = prodsq;
  }
  if (weightname.compare("minsq") == 0){
    h2fun = minsq;
  }
  if (weightname.compare("prod1") == 0){
    h2fun = hprod;
  }
  //check weight function
  if (h2fun == nullptr){
    throw std::invalid_argument("Matching weight function not found");
  }

  //choose manifold object
  manifold<a1type> manobj;
  if (manifoldname.compare("sphere") == 0){
    manobj = Spos;
  }
  if (manifoldname.compare("simplex") == 0){
    manobj = simplex;
  }

  //choose ll function
  a1type (*ll)(const veca1 &, const veca1 &) = nullptr;
  if (llname.compare("dirichlet") == 0){
    ll = ll_dirichlet;
  }
  if (llname.compare("ppi") == 0){
    ll = ll_ppi;
  }
  //check ll function
  if (ll == nullptr){
    throw std::invalid_argument("Matching ll function not found");
  }

  *out = tapesmo(xbetain, n, ll,
                 manobj,
                 h2fun, acut);

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

//calc smo and additions
//use a pointer to an ADFun object to compute the function evaluated  at a location
//' @title The score matching objective calculator.
//' @param u A vector in the simplex.
//' @param betain
//' @return The score matching objective value
//' @export
// [[Rcpp::export]]
double psmo(XPtr< CppAD::ADFun<double> > pfun, svecd u, svecd betain){
  //convert input to an Eigen vectors
  vecd u_e(u.size());
  for (size_t i=0; i<u.size(); i++){
    u_e[i] = u[i];
  }
  vecd beta_e(betain.size());
  for (size_t i=0; i<betain.size(); i++){
    beta_e[i] = betain[i];
  }

  vecd xbetain(u_e.size() + beta_e.size());
  xbetain << u_e, beta_e;
  vecd smo_val(1);
  smo_val = pfun->Forward(0, xbetain);  //treat the XPtr as a regular pointer
  return(smo_val[0]);
}

//calc smo and additions
//use a pointer to an ADFun object to compute the function evaluated  at a location
//' @title The score matching objective calculator.
//' @param u A vector in the simplex.
//' @param betain
//' @return The score matching objective value
//' @export
// [[Rcpp::export]]
svecd psmograd(XPtr< CppAD::ADFun<double> > pfun, svecd u, svecd betain){
  //convert input to an Eigen vectors
  vecd u_e(u.size());
  for (size_t i=0; i<u.size(); i++){
    u_e[i] = u[i];
  }
  vecd beta_e(betain.size());
  for (size_t i=0; i<betain.size(); i++){
    beta_e[i] = betain[i];
  }


  vecd xbetain(u_e.size() + beta_e.size());
  xbetain << u_e, beta_e;
  vecd sc_grad(xbetain.size());
  vecd out_e(beta_e.size());
  svecd out(beta_e.size());
  sc_grad = pfun->Jacobian(xbetain);  //treat the XPtr as a regular pointer
  out_e = sc_grad.block(u_e.size(),0,beta_e.size(),1);

  //convert to std::vector
  for (size_t i = 0; i<u_e.size(); i++){
    out[i] = out_e[i];
  }
  return(out);
}

//for testing
//' @title The ppi likelihood calculation
//' @param u A vector in the simplex.
//' @param beta a vector of parameters
//' @return The loglikelihood value (unnormalised).
//' @export
// [[Rcpp::export]]
double ppill(const svecd &beta,
         const svecd & u){
  //convert input to an Eigen vectors
  vecd u_e(u.size());
  for (size_t i=0; i<u.size(); i++){
    u_e[i] = u[i];
  }
  vecd beta_e(beta.size());
  for (size_t i=0; i<beta.size(); i++){
    beta_e[i] = beta[i];
  }

  double y;
  y = ll_ppi(beta_e, u_e);
  return(y);
}

