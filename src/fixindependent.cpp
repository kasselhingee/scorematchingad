#include "fixindependent.h"

pADFun fixindependent(pADFun & uld, // the unnormalised log density tape
                  veca1 x, //new x to use for taping
                  Eigen::Matrix<int, Eigen::Dynamic, 1> fixedx){ //TRUE (1) values indicate that the corresponding value of x is not a variable (dynamic or independent)
  if (fixedx.size() != uld.size_dyn_ind()){
    Rcpp::stop("fixedx must have the same length as the independent vector of uld");
  }
  if (x.size() != fixedx.size()){
    Rcpp::stop("x and fixedx must have the same length");
  }
  
  //separate fixed and variable x
  veca1 xvar(x.size() - fixedx.sum());
  veca1 xfxd(fixedx.sum());
  size_t idx_var(0);
  size_t idx_fxd(0);
  for (long int i=0;i<x.size();i++){
    if (fixedx[i]){
      xfxd[idx_fxd] = x[i];
      idx_fxd += 1;
    } else {
      xvar[idx_var] = x[i];
      idx_var += 1;
    }
  }

  //convert taped object to execute on a1type (not double)
  CppAD::ADFun<a1type, double> uldhigher;
  uldhigher = (uld.get_ptr())->base2ad();

  //redo tape with some of the dynamic parameters fixed as per fixedx
  veca1 theta = uld.dyntape;
  CppAD::Independent(xvar, theta);


  //combine fixed and variable x
  veca1 xrecom(x.size());
  idx_var = 0;
  for (long int i=0;i<x.size();i++){
    if (fixedx[i]){
      xrecom[i] = x[i];
    } else {
      xrecom[i] = xvar[idx_var];
      idx_var += 1;
    }
  }

  //evaluate tape at new values
  uldhigher.new_dynamic(theta);
  veca1 y(uld.Range());
  y = uldhigher.Forward(0, xrecom);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(x, y);
  tape.check_for_nan(false);

  pADFun out(tape, xvar, theta, uld.name);
  return out;
}

