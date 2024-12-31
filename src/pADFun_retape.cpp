#include "pADFun_retape.h"


pADFun tape_swap(pADFun & pfun){
  //convert taped object to higher order
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = (pfun.get_ptr())->base2ad();

  veca1 y(pfun.Range());
  veca1 newx = pfun.dyntape;
  veca1 newdyn = pfun.xtape;

  //START TAPING
  CppAD::Independent(newx, newdyn);

  pfunhigher.new_dynamic(newx); //before switch the newx is the dynamic parameter vector
  y = pfunhigher.Forward(0, newdyn); //before the switch the newdynparam is the independent value

  //end taping
  CppAD::ADFun<double> tape;
  tape.Dependent(newx, y);
  //out.optimize(); //meant to remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore. But asserts errors.
  tape.check_for_nan(false);

  pADFun out(tape, newx, newdyn, pfun.name);
  return(out);
}


pADFun  tape_Jacobian(pADFun & pfun){
  //convert taped object to higher order, so that the 'base' type of the tape is a1type, so x and dynparam can be passed into Jacobian()
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = (pfun.get_ptr())->base2ad();

  veca1 x = pfun.xtape;
  veca1 dynparam = pfun.dyntape;

  CppAD::Independent(x, dynparam);  //start taping with x as the usual independent parameter and dynparam as the dynamic parameter
  pfunhigher.new_dynamic(dynparam);
  veca1 jac(pfunhigher.Domain() * pfunhigher.Range());
  jac = pfunhigher.Jacobian(x);

  //end taping
  CppAD::ADFun<double> tape;
  tape.Dependent(x, jac);
  //out.optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.  But asserts errors.
  tape.check_for_nan(false);

  pADFun out(tape, x, dynparam, "d(" + pfun.name + ")");
  return(out);
}

pADFun  tape_Hessian(pADFun & pfun){
  if (pfun.Range()>1){
    Rcpp::stop("Taped function 'pfun' must return a vector of length 1. Currently 'pfun' returns a vector of length %i.", pfun.Range());
  }

  //convert taped object to higher order, so that the 'base' type of the tape is a1type, so x and dynparam can be passed into Hessian()
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = (pfun.get_ptr())->base2ad();

  veca1 x = pfun.xtape;
  veca1 dynparam = pfun.dyntape;

  CppAD::Independent(x, dynparam);  //start taping with x as the usual independent parameter and dynparam as the dynamic parameter
  pfunhigher.new_dynamic(dynparam);
  veca1 hess(pfunhigher.Domain() * pfunhigher.Domain());
  hess = pfunhigher.Hessian(x, 0);

  //end taping
  CppAD::ADFun<double> tape;
  tape.Dependent(x, hess);
  //out.optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore. But asserts errors.
  tape.check_for_nan(false);

  pADFun out(tape, x, dynparam, "d^2(" + pfun.name + ")");
  return(out);
}

pADFun  tape_gradoffset(pADFun & pfun){
  if (pfun.Range()>1){
    Rcpp::stop("Taped function 'pfun' must return a vector of length 1. Currently 'pfun' returns a vector of length %i.", pfun.Range());
  }

  //convert taped object to higher order, so that the 'base' type of the tape is a1type, so x and dynparam can be passed into Hessian()
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = (pfun.get_ptr())->base2ad();

  veca1 x = pfun.xtape;
  veca1 dynparam = pfun.dyntape;

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
  CppAD::ADFun<double> tape;
  tape.Dependent(dynparam, gradoffset);
  //out.optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.  But asserts errors.
  tape.check_for_nan(false);

  pADFun out(tape, x, dynparam, "doffset(" + pfun.name + ")");
  return(out);
}

pADFun  tape_logJacdet(pADFun & pfun){
  // domain and range must have equal size for the determinant of the Jacobian to make sense
  if (pfun.Domain() != pfun.Range()){Rcpp::stop("Domain (size %i) and range (size %i) need to be equal for determinant of Jacobian.", pfun.Domain(), pfun.Range());}
  // x and dynparam must have elements of a1type so that taping can proceed
  
  //convert taped object to higher order, so that the 'base' type of the tape is a1type, so x and dynparam can be passed into Jacobian()
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = (pfun.get_ptr())->base2ad();

  veca1 x = pfun.xtape;
  veca1 dynparam = pfun.dyntape;

  CppAD::Independent(x, dynparam);  //start taping with x as the usual independent parameter and dynparam as the dynamic parameter
  pfunhigher.new_dynamic(dynparam);
  mata1 jacmat(pfunhigher.Domain() * pfunhigher.Range(), 1);
  jacmat = pfunhigher.Jacobian(x);
  jacmat.resize(pfunhigher.Domain(), pfunhigher.Range()); //first row is: du1/dz1, du2/dz1, du3/dz1. Second row is du1/dz2, du2/dz2, du3/dz2
  veca1 logdet(1);
  logdet[0] = CppAD::log(CppAD::abs(jacmat.determinant()));

  //end taping
  CppAD::ADFun<double> tape;
  tape.Dependent(x, logdet);
  //out->optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.  But asserts errors.
  tape.check_for_nan(false);

  pADFun out(tape, x, dynparam, "logJdet(" + pfun.name + ")");
  return(out);
}

