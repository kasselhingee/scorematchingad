// Function for fixing some of the dynamic parameters in an existing pADFun

#include "fixdynamic.h"

pADFun fixdynamic(pADFun & pfun, // the unnormalised log density tape
                  veca1 theta, //new theta to use for taping
                  Eigen::Matrix<int, Eigen::Dynamic, 1> isfixed){ //TRUE (1) values indicate that the corresponding value of theta is not a variable (dynamic or independent)
  if (isfixed.size() != pfun.size_dyn_ind()){
    Rcpp::stop("isfixed must have the same length as the dynamic parameter vector of pfun");
  }
  if (theta.size() != isfixed.size()){
    Rcpp::stop("theta and isfixed must have the same length");
  }
  
  //separate fixed and variable theta
  veca1 thetavar(theta.size() - isfixed.sum());
  veca1 thetafxd(isfixed.sum());
  size_t idx_var(0);
  size_t idx_fxd(0);
  for (long int i=0;i<theta.size();i++){
    if (isfixed[i]){
      thetafxd[idx_fxd] = theta[i];
      idx_fxd += 1;
    } else {
      thetavar[idx_var] = theta[i];
      idx_var += 1;
    }
  }

  //convert taped object to execute on a1type (not double)
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = (pfun.get_ptr())->base2ad();

  //redo tape with some of the dynamic parameters fixed as per isfixed
  veca1 x = pfun.xtape;
  CppAD::Independent(x, thetavar);


  //combine fixed and variable theta
  veca1 thetarecom(theta.size());
  idx_var = 0;
  for (long int i=0;i<theta.size();i++){
    if (isfixed[i]){
      thetarecom[i] = theta[i];
    } else {
      thetarecom[i] = thetavar[idx_var];
      idx_var += 1;
    }
  }

  //evaluate tape at new dynamic parameters
  pfunhigher.new_dynamic(thetarecom);
  veca1 y(pfun.Range());
  y = pfunhigher.Forward(0, x);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(x, y);
  tape.check_for_nan(false);

  pADFun out(tape, x, thetavar, pfun.name);
  return out;
}


pADFun reembed(pADFun & uld, 
               transform<a1type> & tran){
  //convert taped object to execute on a1type (not double)
  CppAD::ADFun<a1type, double> uldhigher;
  uldhigher = (uld.get_ptr())->base2ad();

  veca1 z = tran.toM(uld.xtape);
  veca1 theta = uld.dyntape;
  veca1 y(uld.Range());
  veca1 u(uld.xtape.size());


  CppAD::Independent(z, theta);
  u = tran.fromM(z);
  uldhigher.new_dynamic(theta);
  y = uldhigher.Forward(0, u).array() + tran.logdetJfromM(z);

  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(z, y);
  tape.check_for_nan(false);

  pADFun out(tape, z, theta, uld.name+tran.name());
  return out;
}


