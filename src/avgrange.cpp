#include "avgrange.h"

pADFun avgrange(pADFun & pfun){
  //convert taped object to execute on a1type (not double)
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = (pfun.get_ptr())->base2ad();

  veca1 theta = pfun.dyntape;
  veca1 x = pfun.xtape;
  CppAD::Independent(x, theta);

  //evaluate tape
  pfunhigher.new_dynamic(theta);
  veca1 y(pfun.Range());
  y = pfunhigher.Forward(0, x);
  veca1 yavg(1);
  yavg(0) = y.sum()/y.size();
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(x, yavg);
  tape.check_for_nan(false);

  pADFun out(tape, x, theta, "avg_" + pfun.name);
  return out;
}

