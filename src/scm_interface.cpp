#include "scm_interface.hpp"
using namespace Rcpp;

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

veca1 ptoM(XPtr< manifold<a1type> > pman, veca1 u_ad){
  veca1 z_ad(u_ad.size());
  z_ad = pman->toM(u_ad);
  return(z_ad);
}


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


vecd pJacobian(XPtr< CppAD::ADFun<double> > pfun, vecd value, vecd theta){
  //check inputs and tape match
  if (pfun->Domain() != value.size()){stop("Size of input vector %i does not match domain size %i of taped function.", value.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != theta.size()){stop("Size of parameter vector %i does not match parameter size %i of the taped function.", theta.size(), pfun->size_dyn_ind());}

  vecd grad(value.size());
  pfun->new_dynamic(theta);
  grad = pfun->Jacobian(value);  //treat the XPtr as a regular pointer

  return(grad);
}

vecd pForward0(XPtr< CppAD::ADFun<double> > pfun, vecd x, vecd dynparam){
  //check inputs and tape match
  if (pfun->Domain() != x.size()){stop("Size of input vector %i does not match domain size %i of taped function.", x.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != dynparam.size()){stop("Size of parameter vector %i does not match parameter size %i of the taped function.", dynparam.size(), pfun->size_dyn_ind());}

  vecd out(1);
  pfun->new_dynamic(dynparam);
  out = pfun->Forward(0, x);  //treat the XPtr as a regular pointer

  return(out);
}

vecd pHessian(XPtr< CppAD::ADFun<double> > pfun, vecd value, vecd theta){
  //check inputs and tape match
  if (pfun->Domain() != value.size()){stop("Size of input vector %i does not match domain size %i of taped function.", value.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != theta.size()){stop("Size of parameter vector %i does not match parameter size %i of the taped function.", theta.size(), pfun->size_dyn_ind());}

  vecd hess(value.size() * value.size(), 1);
  pfun->new_dynamic(theta);
  hess = pfun->Hessian(value, 0);  //treat the XPtr as a regular pointer
  return(hess);
}


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

std::vector<bool> pParameter(XPtr< CppAD::ADFun<double> > pfun){
  std::vector<bool> isparameter(pfun->Range());
  for (size_t i = 0; i < pfun->Range(); i++){
    isparameter[i] = pfun->Parameter(i);
  }
  return(isparameter);
}
// According to the help, applying Variable(u) to each return value would be false if u depends on the dynamic parameters and does not depend on the independent variable vector.

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

