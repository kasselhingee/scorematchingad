#include "keeprange.h"

pADFun keeprange(pADFun & pfun, // the unnormalised log density tape
                 Eigen::Matrix<int, Eigen::Dynamic, 1> keep){ //indices of range elements to keep. STARTING AT 1 FOR R COMPATABILITY
  if (keep.size() == 0){Rcpp::stop("keep is empty: at least one element of the range must be kept");}
  if (keep.maxCoeff() > pfun.Range()){
    Rcpp::stop("keep has indices larger than range of pfun");
  }
  keep = keep.array() - 1;  //convert to C++ indices

  //convert taped object to execute on a1type (not double)
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = (pfun.get_ptr())->base2ad();

  //redo tape with some of the dynamic parameters fixed as per keep
  veca1 theta = pfun.dyntape;
  veca1 x = pfun.xtape;
  CppAD::Independent(x, theta);

  //evaluate tape at new values
  pfunhigher.new_dynamic(theta);
  veca1 y(pfun.Range());
  y = pfunhigher.Forward(0, x);
  veca1 filtered = y(keep);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(x, filtered);
  tape.check_for_nan(false);

  pADFun out(tape, x, theta, pfun.name);
  return out;
}

