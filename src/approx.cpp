#include "approx.h"
#include <RcppEigen.h>

vecd taylorApprox_currentdynparam(pADFun & pfun,  //a tape with independent values that are points on the manifold (not the parameters)
		  vecd x,
                  vecd centre,
		  const size_t order){
    //In CppAD speak consider the input function X(t) to be
    //X(t) = centre + t*(x - centre). So X(0) = centre, X(1) = x.
    //First derivative of X at 0, is x - centre
    //Higher order derivative are all 0.
    vecd out(pfun.Range());
    vecd diff(x.size());
    out.setZero();
    out += pfun.Forward(0, centre); //zeroth order - constant component of taylor expansion
    if (order >= 1){
      diff = x - centre; // for some reason forward can't take the lazy evaluation version of x - centre direclty. (x - centre).eval() also works
      out += pfun.Forward(1, diff); //now out[0] is evaluation of a linear approximation of f
    }
    if (order >= 2){
      for (size_t i=2; i<=order; i++){
	      diff.setZero();
        out += pfun.Forward(i, diff); //now out[0] is evaluation of a quadratic approximation of f
      }
    }
    return(out);
} 

vecd taylorApprox(pADFun & pfun,  //a tape with independent values that are points on the manifold (not the parameters)
		  vecd x,
                  vecd centre,
                  vecd dynparam,
		  const size_t order){
   pfun.new_dynamic(dynparam);

   return taylorApprox_currentdynparam(pfun, x, centre, order); 
} 

