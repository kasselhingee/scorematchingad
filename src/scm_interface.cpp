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

  pfun->new_dynamic(u_e); //for smo the non-derivative param is u
  vecd smo_val(1);
  smo_val = pfun->Forward(0, beta_e);  //treat the XPtr as a regular pointer
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

  pfun->new_dynamic(u_e); //for smo the non-derivative param is u
  vecd out_e(beta_e.size());
  svecd out(beta_e.size());
  out_e = pfun->Jacobian(beta_e);  //treat the XPtr as a regular pointer

  //convert to std::vector
  for (size_t i = 0; i<out_e.size(); i++){
    out[i] = out_e[i];
  }
  return(out);
}

//' @title Tape of a log-likelihood calculation
//' @param p dimension of measurements
//' @param bd dimension of the parameter vector
//' @param llname name of the likelihood function
//' @return An RCpp::XPtr object pointing to the ADFun
//' @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> > ptapell(
                                      size_t d, //dimensions of measurments
                                      size_t bd, //dimension of parameters
                                      std::string llname){
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

  veca1 utape(d);
  veca1 thetatape(bd);
  utape.setOnes();
  thetatape.setOnes();
  CppAD::ADFun<double> ppitape;

  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  *out = tapell(utape,
                thetatape,
                   ll,
                   simplex::fromM, //transformation from manifold to simplex
                   simplex::logdetJ_fromM); //determinant of Jacobian of the tranformation - for correcting the likelihood function as it is a density

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

//for testing
//' @title Tape a likelihood tape wrt theta
//' @param u A vector in the simplex.
//' @param theta a vector of parameters
//' @return A pointer to an ADFun
//' @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> > ptapell_theta(XPtr< CppAD::ADFun<double> > pfun, svecd u, svecd thetain){
  //convert input to an Eigen vectors
  veca1 u_e(u.size());
  for (size_t i=0; i<u.size(); i++){
    u_e[i] = u[i];
  }
  veca1 theta_e(thetain.size());
  for (size_t i=0; i<thetain.size(); i++){
    theta_e[i] = thetain[i];
  }

  //convert taped object to higher order
  CppAD::ADFun<a1type, double> lltape;
  lltape = pfun->base2ad();

  veca1 y(1);

  //START TAPING
  CppAD::Independent(theta_e, u_e);

  lltape.new_dynamic(theta_e);
  y = lltape.Forward(0, u_e);

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(theta_e, y);
  out->optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.
  out->check_for_nan(false);

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

//for testing
//' @title The Jacobian of recorded function
//' @param pfun Rcpp::XPtr to an ADFun with dynamic parameters
//' @param u A vector in the simplex.
//' @param beta a vector of the dynamic parameters
//' @return The Jacobian of pfun
//' @export
// [[Rcpp::export]]
svecd pJacobian(XPtr< CppAD::ADFun<double> > pfun, svecd u, svecd betain){
  //convert input to an Eigen vectors
  vecd u_e(u.size());
  for (size_t i=0; i<u.size(); i++){
    u_e[i] = u[i];
  }
  vecd beta_e(betain.size());
  for (size_t i=0; i<betain.size(); i++){
    beta_e[i] = betain[i];
  }


  vecd grad(u_e.size());
  svecd out(u_e.size());
  pfun->new_dynamic(beta_e);
  grad = pfun->Jacobian(u_e);  //treat the XPtr as a regular pointer

  //convert to std::vector
  for (size_t i = 0; i<grad.size(); i++){
    out[i] = grad[i];
  }
  return(out);
}

//for testing
//' @title The value of a recorded function
//' @param pfun Rcpp::XPtr to an ADFun with dynamic parameters
//' @param u A vector in the simplex.
//' @param beta a vector of the dynamic parameters
//' @return The value of pfun
//' @export
// [[Rcpp::export]]
double pForward0(XPtr< CppAD::ADFun<double> > pfun, svecd u, svecd betain){
  //convert input to an Eigen vectors
  vecd u_e(u.size());
  for (size_t i=0; i<u.size(); i++){
    u_e[i] = u[i];
  }
  vecd beta_e(betain.size());
  for (size_t i=0; i<betain.size(); i++){
    beta_e[i] = betain[i];
  }


  vecd out_e(1);
  pfun->new_dynamic(beta_e);
  out_e = pfun->Forward(0, u_e);  //treat the XPtr as a regular pointer

  return(out_e[0]);
}

//for testing
//' @title The Hessian of recorded function
//' @param pfun Rcpp::XPtr to an ADFun with dynamic parameters
//' @param u A vector in the simplex.
//' @param beta a vector of the dynamic parameters
//' @return The Hessian of pfun
//' @export
// [[Rcpp::export]]
svecd pHessian(XPtr< CppAD::ADFun<double> > pfun, svecd u, svecd betain){
  //convert input to an Eigen vectors
  vecd u_e(u.size());
  for (size_t i=0; i<u.size(); i++){
    u_e[i] = u[i];
  }
  vecd beta_e(betain.size());
  for (size_t i=0; i<betain.size(); i++){
    beta_e[i] = betain[i];
  }


  vecd hess(u_e.size() * u_e.size(), 1);
  svecd out(hess.size());
  pfun->new_dynamic(beta_e);
  hess = pfun->Hessian(u_e, 0);  //treat the XPtr as a regular pointer

  //convert to std::vector
  for (size_t i = 0; i<hess.size(); i++){
    out[i] = hess[i];
  }
  return(out);
}
