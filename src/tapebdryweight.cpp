# include "tapebdryweight.h"
# include "manifoldtransforms/bdryweights.hpp"

CppAD::ADFun<double> tapeh2(veca1 z,
                            a1type (*h2fun)(const veca1 &, const double &),
                            const double & acut){
  //tape relationship between z and h2
  CppAD::Independent(z);
  // range space vector
  size_t m = 1;               // number of ranges space variables
  veca1 y(m); // vector of ranges space variables
  y[0] = h2fun(z, acut);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(z, y);
  return(tape);
}


pADFun tape_bdryw_inbuilt(std::string name, veca1 x, const double & acut){
  //choose weight function
  a1type (*h2fun)(const veca1 &, const double &) = nullptr;
  if (name.compare("prodsq") == 0){
    h2fun = bdryweight::prodsq;
  }
  if (name.compare("minsq") == 0){
    h2fun = bdryweight::minsq;
  }
  if (name.compare("ones") == 0){
    h2fun = bdryweight::ones;
  }
  //check weight function
  if (h2fun == nullptr){
    throw std::invalid_argument("Matching weight function not found");
  }
  CppAD::ADFun<double> tape;
  CppAD::Independent(x);
  veca1 y(1);
  y(0) = h2fun(x, acut);
  tape.Dependent(x, y);
  tape.check_for_nan(false);
  veca1 empty(0);
  pADFun out(tape, x, empty, name);
  return(out);
}

