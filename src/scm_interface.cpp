//for content that is Rcpp specific

#include "mycpp/scm.cpp"
using namespace Rcpp;

////////////// Create Pointers to Manifold Objects ///////////////
//in R store a pointer to the ADFun object
//' @title Generate manifold object
//' @param manifoldname The name of the manifold to transform to. Either 'sphere' or 'simplex'
//' @return An RCpp::XPtr object pointing to the C++ manifold object
//' @export
// [[Rcpp::export]]
XPtr< manifold<a1type> > pmanifold(std::string manifoldname){
  manifold<a1type> * out = new manifold<a1type>; //returning a pointer
  if (manifoldname.compare("sphere") == 0){
    *out = {
      Spos::toS, Spos::Pmat_S, Spos::dPmat_S,
      Spos::fromS, Spos::logdetJ_fromS,
    };
  }
  if (manifoldname.compare("simplex") == 0){
    *out = {
      simplex::toM, simplex::Pmat_M, simplex::dPmat_M,
      simplex::fromM, simplex::logdetJ_fromM
    };
  }

  XPtr< manifold<a1type> > pout(out, true);
  return(pout);
}

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
XPtr< CppAD::ADFun<double> > ptapesmo(svecd u,
                                      svecd theta,
                                      XPtr< CppAD::ADFun<double> > pll,
                                      XPtr< manifold<a1type> > pman,
                                      std::string weightname,
                                      const double acut,
                                      bool verbose){
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer

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

  //convert svecd to veca1
  veca1 u_ad(u.size());
  for (size_t i=0; i<u.size(); i++){
    u_ad[i] = u[i];
  }
  veca1 theta_ad(theta.size());
  for (size_t i=0; i<theta.size(); i++){
    theta_ad[i] = theta[i];
  }

  *out = tapesmo(u_ad,
                 theta_ad,
                 *pll,
                 *pman,
                 h2fun,
                 acut,
                 verbose);

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

//' @title Tape of a log-likelihood calculation
//' @param p dimension of measurements
//' @param bd dimension of the parameter vector
//' @param llname name of the likelihood function
//' @return An RCpp::XPtr object pointing to the ADFun
//' @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> > ptapell(svecd z, //data measurement on the M manifold
                                     svecd theta,
                                     std::string llname,
                                     XPtr< manifold<a1type> > pman,
                                     svecb fixedtheta,
                                     bool verbose
                                     ){
  Eigen::Matrix<bool, Eigen::Dynamic, 1> fixedtheta_e(fixedtheta.size());
  for (size_t i=0;i<fixedtheta.size();i++){
    fixedtheta_e[i] = fixedtheta[i];
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

  veca1 z_ad(z.size());
  for (size_t i=0; i<z.size(); i++){
    z_ad[i] = z[i];
  }
  veca1 theta_ad(theta.size());
  for (size_t i=0; i<theta.size(); i++){
    theta_ad[i] = theta[i];
  }

  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  *out = tapell(z_ad,
                theta_ad,
                ll,
                pman->fromM, //transformation from manifold to simplex
                pman->logdetJfromM, //determinant of Jacobian of the tranformation - for correcting the likelihood function as it is a density
                fixedtheta_e,
                verbose);

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

//for testing
//' @title Switch Dynamic and pure Independent values
//' @description Convert an ADFun so that the independent values become dynamic parameters
//' and the dynamic parameters become independent values
//' @param newvalue The value (in the sense after the switch has occured) at which to tape the ADFun
//' @param newdynparam The value of the now dynamic parameters at which to tape the ADFun
//' @return A pointer to an ADFun
//' @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> > swapDynamic(XPtr< CppAD::ADFun<double> > pfun, svecd newvalue, svecd newdynparam){
  //convert input to an Eigen vectors
  veca1 value(newvalue.size());
  for (size_t i=0; i<newvalue.size(); i++){
    value[i] = newvalue[i];
  }
  veca1 dynparam(newdynparam.size());
  for (size_t i=0; i<newdynparam.size(); i++){
    dynparam[i] = newdynparam[i];
  }

  //convert taped object to higher order
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = pfun->base2ad();

  veca1 y(1);

  //START TAPING
  CppAD::Independent(value, dynparam);

  pfunhigher.new_dynamic(value); //before switch the value is the dynamic parameter vector
  y = pfunhigher.Forward(0, dynparam); //before the switch the dynparam is the independent value

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(value, y);
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
svecd pJacobian(XPtr< CppAD::ADFun<double> > pfun, svecd value, svecd theta){
  //convert input to an Eigen vectors
  vecd value_e(value.size());
  for (size_t i=0; i<value.size(); i++){
    value_e[i] = value[i];
  }
  vecd theta_e(theta.size());
  for (size_t i=0; i<theta.size(); i++){
    theta_e[i] = theta[i];
  }


  vecd grad(value_e.size());
  svecd out(value_e.size());
  pfun->new_dynamic(theta_e);
  grad = pfun->Jacobian(value_e);  //treat the XPtr as a regular pointer

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
double pForward0(XPtr< CppAD::ADFun<double> > pfun, svecd value, svecd theta){
  //convert input to an Eigen vectors
  vecd value_e(value.size());
  for (size_t i=0; i<value.size(); i++){
    value_e[i] = value[i];
  }
  vecd theta_e(theta.size());
  for (size_t i=0; i<theta.size(); i++){
    theta_e[i] = theta[i];
  }


  vecd out_e(1);
  pfun->new_dynamic(theta_e);
  out_e = pfun->Forward(0, value_e);  //treat the XPtr as a regular pointer

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
svecd pHessian(XPtr< CppAD::ADFun<double> > pfun, svecd value, svecd theta){
  //convert input to an Eigen vectors
  vecd value_e(value.size());
  for (size_t i=0; i<value.size(); i++){
    value_e[i] = value[i];
  }
  vecd theta_e(theta.size());
  for (size_t i=0; i<theta.size(); i++){
    theta_e[i] = theta[i];
  }


  vecd hess(value_e.size() * value_e.size(), 1);
  svecd out(hess.size());
  pfun->new_dynamic(theta_e);
  hess = pfun->Hessian(value_e, 0);  //treat the XPtr as a regular pointer

  //convert to std::vector
  for (size_t i = 0; i<hess.size(); i++){
    out[i] = hess[i];
  }
  return(out);
}
