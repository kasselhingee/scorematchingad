# ifndef SCM_INTERFACE
# define SCM_INTERFACE

//for content that is Rcpp specific
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "mycpp/scm.cpp"
using namespace Rcpp;
#include "mycpp/wrapas.cpp"

////////////// Create Pointers to Manifold Objects ///////////////
//in R store a pointer to the ADFun object
// @title Generate manifold object
//' @param manifoldname The name of the manifold to transform to. Either 'sphere' or 'simplex'
//' @return An RCpp::XPtr object pointing to the C++ manifold object
// @export
// [[Rcpp::export]]
XPtr< manifold<a1type> > pmanifold(std::string manifoldname){
  manifold<a1type> * out;  //returning a pointer
  if (manifoldname.compare("sphere") == 0){
    out = new Spos<a1type>();
  } else if (manifoldname.compare("simplex") == 0){
    out = new simplex<a1type>();
  } else if (manifoldname.compare("Ralr") == 0){
    out = new Ralr<a1type>();
  } else if (manifoldname.compare("Snative") == 0){
    out = new Snative<a1type>();
  } else {
    stop("Manifold not found");
  }

  XPtr< manifold<a1type> > pout(out, true);
  return(pout);
}

// @title Test a manifold object
//' @description A lightweight test of a manifold object.
//' Its main benefit is to force compilation of templated functions for the manifold,
//' and to print results to standard output.
//' @param pman An XPtr to a manifold object. Created by `pmanifold()`
//' @return An integer. 0 if the testable parts pass.
// @export
// [[Rcpp::export]]
int testmanifold(XPtr< manifold<a1type> > pman, veca1 u_ad){
  Rcout << "Starting tests" << std::endl;
  // toM then fromM get back to u
  Rcout << "               Input u was: " << u_ad.transpose() << std::endl;
  veca1 z_ad(u_ad.size());
  z_ad = pman->toM(u_ad);
  Rcout << "                 After toM: " << z_ad.transpose() << std::endl;
  veca1 u2_ad(u_ad.size());
  u2_ad = pman->fromM(z_ad);
  Rcout << "      After toM then fromM: " << u2_ad.transpose() << std::endl;
  if ((u2_ad - u_ad).array().abs().maxCoeff() > 1E-8){
    Rcout << "toM then fromM not passed." << std::endl;
    return(1);
  }

  // Run the other elements
  Rcout << " logdetJ_fromM at toM(u): " << pman->logdetJfromM(z_ad) << std::endl;
  Rcout << " Pmat at toM(u): " << std::endl << pman->Pmatfun(z_ad) << std::endl;
  for (long int d=0; d<u_ad.size(); d++){
    Rcout << " dPmat at toM(u) in dimension " << d <<":" << std::endl << pman->dPmatfun(z_ad, d) << std::endl;
  }
  return(0);
}


// @title Apply to `toM` function of a manifold object
//' @description Apply the `toM` function of a manifold object.
//' @param pman An XPtr to a manifold object. Created by `pmanifold()`.
//' @param u A vector to be transformed to the manifold via `toM`.
//' @return A vector on the manifold.
// @export
// [[Rcpp::export]]
veca1 ptoM(XPtr< manifold<a1type> > pman, veca1 u_ad){
  veca1 z_ad(u_ad.size());
  z_ad = pman->toM(u_ad);
  return(z_ad);
}


//in R store a pointer to the ADFun object
// @title The score matching objective calculator.
//' @param xbetain a concatenated vector of sqrt(x) and beta
//' @param n The dimension of x.
//' @param manifoldname The name of the manifold to transform to
//' @param weightname The name of the weight function to use
//' @param acut The constraint a_c in the weight function
//' @return An RCpp::XPtr object pointing to the ADFun
// @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> > ptapesmo(veca1 u_ad,
                                      veca1 theta_ad,
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
  if (weightname.compare("ones") == 0){
    h2fun = oneweights;
  }
  if (weightname.compare("prod1") == 0){
    h2fun = hprod;
  }
  //check weight function
  if (h2fun == nullptr){
    throw std::invalid_argument("Matching weight function not found");
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

// @title Tape of a log-likelihood calculation
//' @param p dimension of measurements
//' @param bd dimension of the parameter vector
//' @param llname name of the likelihood function
//' @return An RCpp::XPtr object pointing to the ADFun
// @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> > ptapell(veca1 z_ad, //data measurement on the M manifold
                                     veca1 theta_ad,
                                     std::string llname,
                                     XPtr< manifold<a1type> > pman,
                                     Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta,
                                     bool verbose
                                     ){

  //choose ll function
  a1type (*ll)(const veca1 &, const veca1 &) = nullptr;
  if (llname.compare("dirichlet") == 0){
    ll = ll_dirichlet;
  }
  if (llname.compare("ppi") == 0){
    ll = ll_ppi;
  }
  if (llname.compare("vMF") == 0){
    ll = ll_vMF;
  }
  if (llname.compare("Bingham") == 0){
    ll = ll_Bingham;
  }
  if (llname.compare("FB") == 0){
    ll = ll_FB;
  }
  if (llname.compare("Rivest") == 0){
    ll = ll_Rivest;
  }
  //check ll function
  if (ll == nullptr){
    throw std::invalid_argument("Matching ll function not found");
  }


  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  *out = tapell(z_ad,
                theta_ad,
                ll,
                pman.checked_get(),
                fixedtheta,
                verbose);

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

//for testing
// @title Switch Dynamic and pure Independent values
//' @description Convert an ADFun so that the independent values become dynamic parameters
//' and the dynamic parameters become independent values
//' @param newvalue The value (in the sense after the switch has occured) at which to tape the ADFun
//' @param newdynparam The value of the now dynamic parameters at which to tape the ADFun
//' @return A pointer to an ADFun
// @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> > swapDynamic(XPtr< CppAD::ADFun<double> > pfun, veca1 newvalue, veca1 newdynparam){
  //check inputs and tape match
  if (pfun->Domain() != newdynparam.size()){stop("Size of newdynparam must match domain size of taped function.");}
  if (pfun->size_dyn_ind() != newvalue.size()){stop("Size of newvalue must match the parameter size of the taped function.");}



  //convert taped object to higher order
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = pfun->base2ad();

  veca1 y(1);

  //START TAPING
  CppAD::Independent(newvalue, newdynparam);

  pfunhigher.new_dynamic(newvalue); //before switch the newvalue is the dynamic parameter vector
  y = pfunhigher.Forward(0, newdynparam); //before the switch the newdynparam is the independent value

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(newvalue, y);
  out->optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.
  out->check_for_nan(false);

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}


//for testing
// @title The Jacobian of recorded function
//' @param pfun Rcpp::XPtr to an ADFun with dynamic parameters
//' @param u A vector in the simplex.
//' @param beta a vector of the dynamic parameters
//' @return The Jacobian of pfun
// @export
// [[Rcpp::export]]
vecd pJacobian(XPtr< CppAD::ADFun<double> > pfun, vecd value, vecd theta){
  //check inputs and tape match
  if (pfun->Domain() != value.size()){stop("Size of input vector %i does not match domain size %i of taped function.", value.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != theta.size()){stop("Size of parameter vector %i does not match parameter size %i of the taped function.", theta.size(), pfun->size_dyn_ind());}

  vecd grad(value.size());
  pfun->new_dynamic(theta);
  grad = pfun->Jacobian(value);  //treat the XPtr as a regular pointer

  return(grad);
}

// @title The value of a recorded function
//' @param pfun Rcpp::XPtr to an ADFun with dynamic parameters
//' @param value A vector in the domain of the taped function.
//' @param theta a vector of the dynamic parameters
//' @return The value of pfun
// @export
// [[Rcpp::export]]
vecd pForward0(XPtr< CppAD::ADFun<double> > pfun, vecd value, vecd theta){
  //check inputs and tape match
  if (pfun->Domain() != value.size()){stop("Size of input vector %i does not match domain size %i of taped function.", value.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != theta.size()){stop("Size of parameter vector %i does not match parameter size %i of the taped function.", theta.size(), pfun->size_dyn_ind());}

  vecd out(1);
  pfun->new_dynamic(theta);
  out = pfun->Forward(0, value);  //treat the XPtr as a regular pointer

  return(out);
}

// @title The Hessian of recorded function
//' @param pfun Rcpp::XPtr to an ADFun with dynamic parameters
//' @param u A vector in the simplex.
//' @param beta a vector of the dynamic parameters
//' @return The Hessian of pfun
// @export
// [[Rcpp::export]]
vecd pHessian(XPtr< CppAD::ADFun<double> > pfun, vecd value, vecd theta){
  //check inputs and tape match
  if (pfun->Domain() != value.size()){stop("Size of input vector %i does not match domain size %i of taped function.", value.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != theta.size()){stop("Size of parameter vector %i does not match parameter size %i of the taped function.", theta.size(), pfun->size_dyn_ind());}

  vecd hess(value.size() * value.size(), 1);
  pfun->new_dynamic(theta);
  hess = pfun->Hessian(value, 0);  //treat the XPtr as a regular pointer
  return(hess);
}


//' @title The value of a recorded function approximated by Taylor expansion
//' @param pfun Rcpp::XPtr to an ADFun tape a tape with independent values that are the points to be differentiated with
//' @param u A vector in the domain of the taped function.
//' @param centre A vector in the domain of the taped function to approximate the value at `u` from.
//' @param dynparam a vector of the dynamic parameters
//' @param order The order of Taylor expansion to use.
//' @description Approximates the value of a `CppAD` tape at `u` using a Taylor approximation at `centre`. The dynamic parameters of the tape are set by `dynparam`.
//' @return The approximate value of pfun
//' @export
// [[Rcpp::export]]
vecd pTaylorApprox(XPtr< CppAD::ADFun<double> > pfun,
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

// @title The approximate value of the gradient (wrt space 1) of recorded function
//' @param pfun Rcpp::XPtr to an ADFun tape a tape with dynamic parameters and independent parameters
//' @param value A vector in the domain of the taped function.
//' @param thetacentre A vector in the space of the dynamic parameters of the recorded function
//' this vector forms the centre of the Taylor approximation
//' @param theta a vector of the dynamic parameters
//' @param order The order of Taylor expansion to use.
//' @description Taylor expansion in the `theta` dimensions, to approximate the gradient wrt the `value` dimensions.
//' @return The approximate value of the gradient, with respect to theta, of pfun
// @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> >  pTapeJacobianSwap(XPtr< CppAD::ADFun<double> > pfun,
                    veca1 value, veca1 theta){

  //convert taped object to higher order
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = pfun->base2ad();

  //first tape the Jacobian of pfun but with the theta becoming the independent variables
  CppAD::Independent(theta, value);  //for this tape, theta must be altered using new_dynamic
  pfunhigher.new_dynamic(theta);
  veca1 grad(value.size());
  grad = pfunhigher.Jacobian(value);

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(theta, grad);
  out->optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.
  out->check_for_nan(false);

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

// @title The approximate value of the gradient (wrt space 1) of recorded function
//' @param pfun Rcpp::XPtr to an ADFun tape a tape with dynamic parameters and independent parameters
//' @param value A vector in the domain of the taped function.
//' @param thetacentre A vector in the space of the dynamic parameters of the recorded function
//' this vector forms the centre of the Taylor approximation
//' @param theta a vector of the dynamic parameters
//' @param order The order of Taylor expansion to use.
//' @description Taylor expansion in the `theta` dimensions, to approximate the gradient wrt the `value` dimensions.
//' @return The approximate value of the gradient, with respect to theta, of pfun
// @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> >  pTapeHessianSwap(XPtr< CppAD::ADFun<double> > pfun,
                                                veca1 value, veca1 theta){
  //convert taped object to higher order
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = pfun->base2ad();

  //first tape the Jacobian of pfun but with the theta becoming the independent variables
  CppAD::Independent(theta, value);  //for this tape, theta must be altered using new_dynamic
  pfunhigher.new_dynamic(theta);
  veca1 hess(value.size() * value.size());
  hess = pfunhigher.Hessian(value, 0);

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(theta, hess);
  out->optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.
  out->check_for_nan(false);

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

//' @title Tape the Jacobian of CppAD Tape
//' @param pfun Rcpp::XPtr to an ADFun tape a tape with dynamic parameters and independent parameters
//' @param x A vector in the domain of the taped function.
//' this vector forms the centre of the Taylor approximation
//' @param dynparam a vector of the dynamic parameters
//' @description Creates a tape of the Jacobian of function taped by CppAD.
//' When the function returns a real value (as is the case for densities and the score matching objective) the Jacobian is equivalent to the gradient.
//' The `x` vector is used as the value to conduct the taping.
//' @details
//' When the returned tape is evaluated (via say [`pForward0()`], the resultant vector contains the Jacobian in long format (see [https://coin-or.github.io/CppAD/doc/jacobian.htm]).
//' Suppose the function represented by `pfun` maps from \eqn{n}-dimensional space to \eqn{m}-dimensional space, then
//' the first \eqn{n} elements of vector is the gradient of the first component of function output.
//' The next \eqn{n} elements of the vector is the gradient of the second component of the function output.
//' The Jacobian as a matrix, could then be obtained by [`as.matrix()`] with `byrow = TRUE` and `ncol = n`.
//' @return A `Rcpp::XPtr` to a CppAD::ADFun object.
//' @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> >  pTapeJacobian(XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam){
  // x and dynparam must have elements of a1type so that taping can proceed


  //convert taped object to higher order, so that the 'base' type of the tape is a1type, so x and dynparam can be passed into Jacobian()
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = pfun->base2ad();

  CppAD::Independent(x, dynparam);  //start taping with x as the usual independent parameter and dynparam as the dynamic parameter
  pfunhigher.new_dynamic(dynparam);
  veca1 jac(pfunhigher->Domain() * pfunhigher->Range());
  jac = pfunhigher.Jacobian(x);

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(jac, dynparam);
  out->optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.
  out->check_for_nan(false);

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

# endif
