#include "fixindependent.h"

pADFun fixindependent(pADFun & pfun, // the unnormalised log density tape
                  veca1 x, //new x to use for taping
                  Eigen::Matrix<int, Eigen::Dynamic, 1> isfixed){ //TRUE (1) values indicate that the corresponding value of x is not a variable (dynamic or independent)
  if (isfixed.size() != pfun.Domain()){
    Rcpp::stop("isfixed must have the same length as the domain of pfun");
  }
  if (x.size() != isfixed.size()){
    Rcpp::stop("x and isfixed must have the same length");
  }
  
  //separate fixed and variable x
  veca1 xvar(x.size() - isfixed.sum());
  veca1 xfxd(isfixed.sum());
  size_t idx_var(0);
  size_t idx_fxd(0);
  for (long int i=0;i<x.size();i++){
    if (isfixed[i]){
      xfxd[idx_fxd] = x[i];
      idx_fxd += 1;
    } else {
      xvar[idx_var] = x[i];
      idx_var += 1;
    }
  }

  //convert taped object to execute on a1type (not double)
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = (pfun.get_ptr())->base2ad();

  //redo tape with some of the dynamic parameters fixed as per isfixed
  veca1 theta = pfun.dyntape;
  CppAD::Independent(xvar, theta);


  //combine fixed and variable x
  veca1 xrecom(x.size());
  idx_var = 0;
  for (long int i=0;i<x.size();i++){
    if (isfixed[i]){
      xrecom[i] = x[i];
    } else {
      xrecom[i] = xvar[idx_var];
      idx_var += 1;
    }
  }

  //evaluate tape at new values
  pfunhigher.new_dynamic(theta);
  veca1 y(pfun.Range());
  y = pfunhigher.Forward(0, xrecom);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(xvar, y);
  tape.check_for_nan(false);

  pADFun out(tape, xvar, theta, pfun.name);
  return out;
}

