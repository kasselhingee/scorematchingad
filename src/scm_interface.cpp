# ifndef SCM_INTERFACE
# define SCM_INTERFACE

//for content that is Rcpp specific
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
#include "tapell.hpp"
#include "tapesmo.hpp"
#include "mycpp/mantrans.hpp"
#include "mycpp/divweights.hpp"
#include "mycpp/likelihoods.hpp"
#include "mycpp/approx.hpp"
#include "mycpp/wrapas.cpp"

////////////// Create Pointers to Manifold Objects ///////////////
//in R store a pointer to the ADFun object
//' @noRd
//' @title Generate manifold with transformation object
//' @param manifoldname The name of the manifold to transform to. Either 'sphere' or 'simplex'
//' @return An RCpp::XPtr object pointing to the C++ manifold object
//' @details
//'  + "sphere" for square-root transformation from the simplex to the positive orthant of the sphere
//'  + "simplex" for the simplex without any transformation.
//'  + "Ralr" for the additive log-ratio transformation from the simplex to Euclidean space, using the final component of vectors in the denominator of the ratio.
//'  + "Snative" for the sphere without any transformation
//' @export
// [[Rcpp::export]]
XPtr< manifold<a1type> > pmanifold(std::string manifoldname){
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
    stop("Manifold not found");
  }

  XPtr< manifold<a1type> > pout(out, true);
  return(pout);
}

//' @noRd
//' @title Test a manifold object
//' @description A lightweight test of a manifold object.
//' Its main benefit is to force compilation of templated functions for the manifold,
//' and to print results to standard output.
//' @param pman An XPtr to a manifold object. Created by `pmanifold()`
//' @return An integer. 0 if the testable parts pass.
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

//' @noRd
//' @title Apply to `toM` function of a manifold object
//' @description Apply the `toM` function of a manifold object.
//' @param pman An XPtr to a manifold object. Created by `pmanifold()`.
//' @param u A vector to be transformed to the manifold via `toM`.
//' @return A vector on the manifold.
// [[Rcpp::export]]
veca1 ptoM(XPtr< manifold<a1type> > pman, veca1 u_ad){
  veca1 z_ad(u_ad.size());
  z_ad = pman->toM(u_ad);
  return(z_ad);
}


//in R store a pointer to the ADFun object
//' @noRd
//' @title The score matching objective calculator.
//' @param xbetain a concatenated vector of sqrt(x) and beta
//' @param n The dimension of x.
//' @param manifoldname The name of the manifold to transform to
//' @param weightname The name of the weight function to use
//' @param acut The constraint a_c in the weight function
//' @return An RCpp::XPtr object pointing to the ADFun
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
    h2fun = divweight::prodsq;
  }
  if (weightname.compare("minsq") == 0){
    h2fun = divweight::minsq;
  }
  if (weightname.compare("ones") == 0){
    h2fun = divweight::oneweights;
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

//' @noRd
//' @title Tape of a log-likelihood calculation
//' @param p dimension of measurements
//' @param bd dimension of the parameter vector
//' @param llname name of the likelihood function
//' @return An RCpp::XPtr object pointing to the ADFun
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
  if (llname.compare("Rivest") == 0){
    ll = ll::ll_Rivest;
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

//' @title Switch Dynamic and Independent Values of a Tape
//' @description Convert an ADFun so that the independent values become dynamic parameters
//' and the dynamic parameters become independent values
//' @param newvalue The independent value (in the sense after the switch has occurred) at which to tape the ADFun
//' @param newdynparam The value of the dynamic parameters (after the switch) at which to tape the ADFun
//' @return A pointer to an ADFun
//' @export
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


//' @title Evaluate the Jacobian of a tape
//' @param pfun Rcpp::XPtr to an ADFun
//' @param x A vector in the domain of the taped function.
//' @param dynparam a vector of the dynamic parameters. If `pfun` has no dynamic parameters then set `dynparam = vector(mode = "numeric")`.
//' @return The Jacobian of pfun
//' @export
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

//' @title Evaluate a CppAD tape
//' @param pfun Rcpp::XPtr to an ADFun with dynamic parameters
//' @param x A vector in the domain of the taped function.
//' @param dynparam a vector of the dynamic parameters.
//' @return The value of `pfun` evaluated at `x` with parameters `dynparam`.
//' @export
// [[Rcpp::export]]
vecd pForward0(XPtr< CppAD::ADFun<double> > pfun, vecd x, vecd dynparam){
  //check inputs and tape match
  if (pfun->Domain() != x.size()){stop("Size of input vector %i does not match domain size %i of taped function.", x.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != dynparam.size()){stop("Size of parameter vector %i does not match parameter size %i of the taped function.", dynparam.size(), pfun->size_dyn_ind());}

  vecd out(1);
  pfun->new_dynamic(dynparam);
  out = pfun->Forward(0, x);  //treat the XPtr as a regular pointer

  return(out);
}

//' @noRd
//' @title OBSOLETE: The Hessian of recorded function. Used only in smobj.R
//' @param pfun Rcpp::XPtr to an ADFun with dynamic parameters
//' @param u A vector in the simplex.
//' @param beta a vector of the dynamic parameters
//' @return The Hessian of pfun
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

//' @title Tape the Jacobian of CppAD Tape
//' @param pfun Rcpp::XPtr to an ADFun tape a tape with dynamic parameters and independent parameters
//' @param x A vector in the domain of the taped function.
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
  //check inputs and tape match
  if (pfun->Domain() != x.size()){stop("Size of input vector %i does not match domain size %i of taped function.", x.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != dynparam.size()){stop("Size of parameter vector %i does not match parameter size %i of the taped function.", dynparam.size(), pfun->size_dyn_ind());}



  //convert taped object to higher order, so that the 'base' type of the tape is a1type, so x and dynparam can be passed into Jacobian()
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = pfun->base2ad();

  CppAD::Independent(x, dynparam);  //start taping with x as the usual independent parameter and dynparam as the dynamic parameter
  pfunhigher.new_dynamic(dynparam);
  veca1 jac(pfunhigher.Domain() * pfunhigher.Range());
  jac = pfunhigher.Jacobian(x);

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(x, jac);
  out->optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.
  out->check_for_nan(false);

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

//' @title Tape the Hessian of a CppAD Tape
//' @inheritParams pTapeJacobian
//' @description Creates a tape of the Hessian of a function taped by CppAD.
//' The taped function represented by `pfun` must be scalar-valued (i.e. a vector of length 1).
//' The `x` vector and `dynparam` are used as the values to conduct the taping.
//' @details
//' When the returned tape is evaluated (via say [`pForward0()`], the resultant vector contains the Hessian in long format (see [https://coin-or.github.io/CppAD/doc/hessian.htm]).
//' Suppose the function represented by `pfun` maps from \eqn{n}-dimensional space to \eqn{1}-dimensional space, then
//' the first \eqn{n} elements of the vector is the gradient of the partial derivative with respect to the first dimension of the function's domain.
//' The next \eqn{n} elements of the vector is the gradient of the partial derivative of the second dimension of the function's domain.
//' The Hessian as a matrix, can be obtained by using [`as.matrix()`] with `ncol = n`.
//' @return A `Rcpp::XPtr` to a CppAD::ADFun object.
//' @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> >  pTapeHessian(XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam){
  // x and dynparam must have elements of a1type so that taping can proceed
  //check inputs and tape match
  if (pfun->Domain() != x.size()){stop("Size of input vector %i does not match domain size %i of taped function.", x.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != dynparam.size()){stop("Size of parameter vector %i does not match parameter size %i of the taped function.", dynparam.size(), pfun->size_dyn_ind());}


  if (pfun->Range()>1){
    stop("Taped function 'pfun' must return a vector of length 1. Currently 'pfun' returns a vector of length %i.", pfun->Range());
  }

  //convert taped object to higher order, so that the 'base' type of the tape is a1type, so x and dynparam can be passed into Hessian()
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = pfun->base2ad();

  CppAD::Independent(x, dynparam);  //start taping with x as the usual independent parameter and dynparam as the dynamic parameter
  pfunhigher.new_dynamic(dynparam);
  veca1 hess(pfunhigher.Domain() * pfunhigher.Domain());
  hess = pfunhigher.Hessian(x, 0);

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(x, hess);
  out->optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.
  out->check_for_nan(false);

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

//' @title Indicate Constant Components of Range
//' @description Use `CppAD`'s `Parameter()` function for `ADFun` objects to see if the returned values of a tape are constant with respect to the independent values.
//' @param pfun A CppAD tape.
//' @param dynparam A set of dynamic parameters for `pfun`.
//' @return A vector logical values. `TRUE` indicates that element of the tape result is constant.
//' @details The `CppAD` function `Parameter(i)` [https://coin-or.github.io/CppAD/doc/fun_property.htm] returns `TRUE` when the `i`th component of the range does not depend on the independent value
//' (the `i`th component may still depend on the value of the dynamic parameters (see 'Dynamic' in [https://coin-or.github.io/CppAD/doc/glossary.htm#Parameter]) ).
//' @export
// [[Rcpp::export]]
std::vector<bool> pParameter(XPtr< CppAD::ADFun<double> > pfun){
  std::vector<bool> isparameter(pfun->Range());
  for (size_t i = 0; i < pfun->Range(); i++){
    isparameter[i] = pfun->Parameter(i);
  }
  return(isparameter);
}
// According to the help, applying Variable(u) to each return value would be false if u depends on the dynamic parameters and does not depend on the independent variable vector.

//' @title Tape the Gradient Offset of a Quadratic CppAD Tape
//' @inheritParams pTapeJacobian
//' @description A quadratic function can be written as
//' \deqn{f(x;\theta) = \frac{1}{2} x^T W(\theta) x + b(\theta)^Tx + c.}
//' The function `pTapeGradOffset` creates a tape of \eqn{b(\theta)} where \eqn{\theta} is the independent variable.
//' @param pfun A quadratic CppAD Tape. Test for quadratic form using [`testquadratictape()`].
//' @details
//' The gradient of \eqn{f(x; \theta)} with respect to \eqn{\theta} is
//' \deqn{\Delta f(x; \theta) = \frac{1}{2}(W(\theta) + W(\theta)^T)\theta + b(\theta),}
//' and the Hessian is 
//' \deqn{H f(x; \theta) = \frac{1}{2}(W(\theta) + W(\theta)^T),}
//' which does not depend on \eqn{x}.
//' The gradient of the function can be rewritten as
//' \deqn{\Delta f(x;\theta) = H f(x; \theta) x + b(\theta)^T x.}
//' The tape calculates \eqn{b(\theta)} as
//'  \deqn{b(\theta) = \Delta f(x;\theta) - H f(x; \theta) x},
//' which does not depend on \eqn{x}.
//' In `pTapeGradOffset()` the `x` provided as an argument is used as the template for calculating \eqn{b(\theta)}.
//' The `x` vector and `dynparam` are used as the values to conduct the taping.
//' @return A `Rcpp::XPtr` to a CppAD::ADFun object. The independent argument to the function are the dynamic parameters of `pfun`.
//' @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> >  pTapeGradOffset(XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam){
  // x and dynparam must have elements of a1type so that taping can proceed
  //check inputs and tape match
  if (pfun->Domain() != x.size()){stop("Size of input vector %i does not match domain size %i of taped function.", x.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != dynparam.size()){stop("Size of parameter vector %i does not match parameter size %i of the taped function.", dynparam.size(), pfun->size_dyn_ind());}


  if (pfun->Range()>1){
    stop("Taped function 'pfun' must return a vector of length 1. Currently 'pfun' returns a vector of length %i.", pfun->Range());
  }

  //convert taped object to higher order, so that the 'base' type of the tape is a1type, so x and dynparam can be passed into Hessian()
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = pfun->base2ad();

  CppAD::Independent(dynparam);  
  pfunhigher.new_dynamic(dynparam);
  veca1 jac(pfunhigher.Domain());
  jac = pfunhigher.Jacobian(x);
  mata1 hess(pfunhigher.Domain() * pfunhigher.Domain(), 1);
  hess = pfunhigher.Hessian(x, 0);
  //arrange hess into a matrix
  hess.resize(pfunhigher.Domain(),pfunhigher.Domain());

  veca1 gradoffset(pfunhigher.Domain());
  gradoffset = jac - (hess * x);

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(dynparam, gradoffset);
  out->optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.
  out->check_for_nan(false);

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

# endif
